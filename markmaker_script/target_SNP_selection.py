import os
import subprocess as sbp
from datetime import datetime

class target_SNP_selection(object):
    def __init__(self,abam,bbam,cbam,aname,bname,csim,reference,position,output_dir,min_depth,max_depth,bhetero,bsim,minMQ,minBQ):

        self.abam,self.bbam,self.cbam,self.aname,self.bname,self.csim = abam,bbam,cbam,aname,bname,csim
        self.reference,self.position,self.output_dir,self.bhetero,self.bsim = reference,position,output_dir,bhetero,bsim
        self.min_depth,self.max_depth,self.minMQ,self.minBQ = min_depth,max_depth,minMQ,minBQ

        position_array=self.position.split(':')
        self.chr=position_array[0]
        self.s_posi=position_array[1]
        self.e_posi=position_array[2]

    def preparation(self):
        # ディレクトリ作成
        os.makedirs('{0}'.format(self.output_dir), exist_ok=True)
        os.makedirs('{0}/log'.format(self.output_dir), exist_ok=True)
        os.makedirs('{0}/target_SNP_selection'.format(self.output_dir), exist_ok=True)
        # setting情報を保存
        output_set_file = open('{0}/log/setting.log'.format(self.output_dir), "w")
        output_set_file.write('A name={0}\nB name={1}\nA bam={2}\nB bam={3}\nC bam={4}\n'.format(self.aname,self.bname,self.abam,self.bbam,self.cbam))
        output_set_file.write('C sim file={0}\nReference fasta={1}\nTarget region={2}:{3}-{4}\n'.format(self.csim,self.reference,self.chr,self.s_posi,self.e_posi))
        output_set_file.write('Min depth={0}\nMax depth={1}\nMin MQ={2}\nMin BQ={3}\n'.format(self.min_depth,self.max_depth,self.minMQ,self.minBQ))
        output_set_file.write('B select hetero={0}\nB sim file={1}\n'.format(self.bhetero,self.bsim))
        output_set_file.close()

    def callSNP(self):
        def stat_format_change(ref_base,stat):
            stat=stat.upper()
            stat=stat.replace(",",ref_base)
            stat=stat.replace(".",ref_base)
            return stat

        def depth_checker(depth):
            depth_error="."
            if depth==0:
                depth_error="0"
            elif self.min_depth>depth:
                depth_error="L"
            elif self.max_depth<depth:
                depth_error="H"
            return depth_error

        def indel_checker(stat):
            variant_error="."
            if "+" in stat:
                variant_error="+"
            elif "-" in stat:
                variant_error="-"
            return variant_error

        def allele_counter(stat):
            allele_count=0
            for i in ["A","T","G","C"]:
                if stat.count(i)>0:
                    allele_count=allele_count+1
            return allele_count

        def calculation_index(base,stat):
            base_count1=0
            base_count2=0
            SNP_base=base
            for i in ["A","T","G","C"]:
                if base==i:
                    base_count1=stat.count(i)
                elif stat.count(i)>0:
                    SNP_base=i
                    base_count2=stat.count(i)
            index=0
            if base_count1+base_count2>0:
                index=int(base_count2*1000/(base_count1+base_count2))/1000
            return SNP_base,index

        # outputするファイルを用意
        output_log_file = open('{0}/log/samtools.log'.format(self.output_dir), "w")
        output_varinat_file = open('{0}/target_SNP_selection/variant.txt'.format(self.output_dir), "w")
        output_depth_file = open('{0}/target_SNP_selection/error_depth.txt'.format(self.output_dir), "w")
        output_SNP_file = open('{0}/target_SNP_selection/SNP.txt'.format(self.output_dir), "w")
        output_seqA_file = open('{0}/target_SNP_selection/{1}.{2}.{3}-{4}.seq'.format(self.output_dir,self.aname,self.chr,self.s_posi,self.e_posi), "w")
        output_seqB_file = open('{0}/target_SNP_selection/{1}.{2}.{3}-{4}.seq'.format(self.output_dir,self.bname,self.chr,self.s_posi,self.e_posi), "w")
        output_seqN_file = open('{0}/target_SNP_selection/shared.{2}.{3}-{4}.seq'.format(self.output_dir,self.bname,self.chr,self.s_posi,self.e_posi), "w")

        # bsimを辞書に入れる
        b_dict_bottom={}
        b_dict_top={}
        if self.bhetero == "yes":
            with open(self.bsim) as f:
                for line in f:
                    line=line.strip()
                    array=line.split()
                    depth=int(array[0])
                    bottom=float(array[1])
                    top=float(array[2])
                    b_dict_bottom[depth]=bottom
                    b_dict_top[depth]=top

        # csimを辞書に入れる
        c_dict_bottom={}
        c_dict_top={}
        if self.cbam != "":
            with open(self.csim) as f:
                for line in f:
                    line=line.strip()
                    array=line.split()
                    depth=int(array[0])
                    bottom=float(array[1])
                    top=float(array[2])
                    c_dict_bottom[depth]=bottom
                    c_dict_top[depth]=top
        
        # bcftoolsでvariant call
        cmd=''
        if self.cbam=='':
            cmd = 'samtools mpileup --no-output-ends -a -r {0}:{1}-{2} -q {3} -Q {4} -f {5} {6} {7}'.format(self.chr,self.s_posi,self.e_posi,self.minMQ,self.minBQ,self.reference,self.abam,self.bbam)
        else:
            cmd = 'samtools mpileup --no-output-ends -a -r {0}:{1}-{2} -q {3} -Q {4} -f {5} {6} {7} {8}'.format(self.chr,self.s_posi,self.e_posi,self.minMQ,self.minBQ,self.reference,self.abam,self.bbam,self.cbam)
        
        cmd_processsbp=sbp.Popen(cmd,stdout=sbp.PIPE,stderr=output_log_file,shell=True)

        for vcf in cmd_processsbp.stdout:
            line = vcf.decode().strip()
            array=line.split()
            chromosome=array[0]
            posi=int(array[1])
            ref_base=array[2].upper()
            a_DP=int(array[3])
            a_stat=array[4]
            b_DP=int(array[6])
            b_stat=array[7]

            share_base=ref_base

            a_depth_error=depth_checker(a_DP)
            b_depth_error=depth_checker(b_DP)
            if a_depth_error != "." or b_depth_error != ".":
                output_depth_file.write("{0}\t{1}\t{2}\t{3}\n".format(chromosome,posi,a_depth_error,b_depth_error))
                share_base="N"

            a_stat=stat_format_change(ref_base,a_stat)
            a_variant=indel_checker(a_stat)
            a_base=ref_base
            a_index=0
            if a_variant=="." and allele_counter(a_stat)>2:
                a_variant="m"

            b_stat=stat_format_change(ref_base,b_stat)
            b_variant=indel_checker(b_stat)
            b_base=ref_base
            b_index=0
            if b_variant=="." and allele_counter(b_stat)>2:
                b_variant="m"

            if a_variant != "." or b_variant != ".":
                output_varinat_file.write("{0}\t{1}\t{2}\t{3}\n".format(chromosome,posi,a_variant,b_variant))
                share_base="N"
            else:
                merge_stat="{0}{1}".format(a_stat,b_stat)
                if allele_counter(merge_stat)>2:
                    output_varinat_file.write("{0}\t{1}\tM\tM\n".format(chromosome,posi))
                    share_base="N"
                else:
                    a_base,a_index=calculation_index(ref_base,a_stat)
                    if a_index != 0 and a_index != 1:
                        share_base="N"
                        b_base,b_index=calculation_index(ref_base,b_stat)
                        if b_index != 0 and b_index != 1:
                            b_variant="v"
                        output_varinat_file.write("{0}\t{1}\tv\t{2}\n".format(chromosome,posi,b_variant))
                    else:
                        b_base,b_index=calculation_index(a_base,b_stat)
                        if a_index == 1:
                            a_index=0

                        if b_index==1:
                            share_base="N"
                            output_varinat_file.write("{0}\t{1}\ts\ts\n".format(chromosome,posi))
                        elif b_index != 0:
                            share_base="N"
                            output_varinat_file.write("{0}\t{1}\t.\tv\n".format(chromosome,posi))

                        if a_depth_error == "." and b_depth_error == ".":
                            if self.bhetero == "yes":
                                if b_dict_bottom[b_DP]<=b_index and b_dict_top[b_DP]>=b_index:
                                    if self.cbam != "":
                                        c_DP=int(array[9])
                                        c_stat=array[10]
                                        c_stat=stat_format_change(ref_base,c_stat)
                                        if depth_checker(c_DP)=="." and indel_checker(c_stat)==".":
                                                merge_stat="{0}{1}{2}".format(a_stat,b_stat,c_stat)
                                                if allele_counter(merge_stat)==2:
                                                    c_base,c_index=calculation_index(a_base,c_stat)
                                                    if c_dict_bottom[c_DP]<=c_index and c_dict_top[c_DP]>=c_index:
                                                        output_SNP_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chromosome,posi,a_base,0,b_base,b_index))
                                    else:
                                        output_SNP_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chromosome,posi,a_base,0,b_base,b_index))
                            elif b_index==1:
                                if self.cbam != "":
                                    c_DP=int(array[9])
                                    c_stat=array[10]
                                    c_stat=stat_format_change(ref_base,c_stat)
                                    if depth_checker(c_DP)=="." and indel_checker(c_stat)==".":
                                            merge_stat="{0}{1}{2}".format(a_stat,b_stat,c_stat)
                                            if allele_counter(merge_stat)==2:
                                                c_base,c_index=calculation_index(a_base,c_stat)
                                                if c_dict_bottom[c_DP]<=c_index and c_dict_top[c_DP]>=c_index:
                                                    output_SNP_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chromosome,posi,a_base,0,b_base,b_index))
                                else:
                                    output_SNP_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chromosome,posi,a_base,0,b_base,b_index))

            output_seqA_file.write(a_base)
            output_seqB_file.write(b_base)
            output_seqN_file.write(share_base)
            

    def run(self):
        print('[MarkMaker:{}] Start target_SNP_selection'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.preparation()
        print('[MarkMaker:{}] Call SNP'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.callSNP()

        output = open('{0}/log/target_SNP_selection.log'.format(self.output_dir), "w")
        count = 0
        with open('{0}/target_SNP_selection/SNP.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                count += 1
        output.write("{0} SNPs were selected".format(count))

        print('[MarkMaker:{}] Finish target_SNP_selection'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        
