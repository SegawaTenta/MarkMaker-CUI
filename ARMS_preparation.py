import os
import primer3
from datetime import datetime

class ARMS_preparation(object):
    def __init__(self,output_dir,recipe,thread):
        self.output_dir,self.recipe,self.thread=output_dir,recipe,thread
        
    def preparation(self):
        # setting情報から設定を入手
        self.aname=""
        self.bname=""
        self.chr=""
        self.s_posi=0
        self.e_posi=0
        with open('{0}/log/setting.log'.format(self.output_dir)) as f:
            for line_count,line in enumerate(f):
                line=line.strip()
                array=line.split("=")
                info=array[1]
                if line_count==0:
                    self.aname=info
                elif line_count==1:
                    self.bname=info
                elif line_count==7:
                    info2=info.split(":")
                    self.chr=info2[0]
                    info3=info2[1].split("-")
                    self.s_posi=int(info3[0])
                    self.e_posi=int(info3[1])
        # recipeから設定を入手
        self.opt_GC=50
        self.min_GC=40
        self.max_GC=60
        self.opt_ln=20
        self.min_ln=17
        self.max_ln=25
        self.opt_tm=60
        self.min_tm=62.5
        self.max_tm=65

        self.max_poly_x=4
        self.max_self_any=8
        self.max_self_end=3

        with open(self.recipe) as f:
            for line in f:
                line=line.strip()
                array=line.split("=")
                set_name=array[0]
                info=float(array[1])
                if "PRIMER_OPT_SIZE"==set_name:
                    self.opt_ln=int(info)
                elif "PRIMER_MIN_SIZE"==set_name:
                    self.min_ln=int(info)
                elif "PRIMER_MAX_SIZE"==set_name:
                    self.max_ln=int(info)
                elif "PRIMER_OPT_TM"==set_name:
                    self.opt_tm=info
                elif "PRIMER_MIN_TM"==set_name:
                    self.min_tm=info
                elif "PRIMER_MAX_TM"==set_name:
                    self.max_tm=info
                elif "PRIMER_OPT_GC"==set_name:
                    self.opt_GC=info
                elif "PRIMER_MIN_GC"==set_name:
                    self.min_GC=info
                elif "PRIMER_MAX_GC"==set_name:
                    self.max_GC=info
                elif "PRIMER_MAX_POLY_X"==set_name:
                    self.max_poly_x=info
                elif "PRIMER_MAX_SELF_ANY"==set_name:
                    self.max_self_any=info
                elif "PRIMER_MAX_SELF_END"==set_name:
                    self.max_self_end=info

        os.makedirs('{0}/ARMS_preparation'.format(self.output_dir), exist_ok=True)

        self.dict_variant_error={}
        with open('{0}/target_SNP_selection/variant.txt'.format(self.output_dir)) as variant_txt:
            for line in variant_txt:
                line=line.strip()
                array=line.split("\t")
                self.dict_variant_error[array[1]]="{0}:{1}".format(array[2],array[3])
        self.dict_depth_error={}
        with open('{0}/target_SNP_selection/error_depth.txt'.format(self.output_dir)) as variant_txt:
            for line in variant_txt:
                line=line.strip()
                array=line.split("\t")
                self.dict_depth_error[array[1]]="{0}:{1}".format(array[2],array[3])

        self.len_list=[self.opt_ln]
        diff_ln=self.max_ln-self.opt_ln
        if self.opt_ln-self.max_ln>diff_ln:
            diff_ln=self.opt_ln-self.max_ln
        for i in range(1,diff_ln+1):
            if self.opt_ln+i<=self.max_ln:
                self.len_list.append(self.opt_ln+i)
            if self.opt_ln-i>=self.min_ln:
                self.len_list.append(self.opt_ln-i)

    def make_primer(self):
        output_test_primers = open('{0}/ARMS_preparation/made_primers.txt'.format(self.output_dir), "w")

        def run_primer3(sequence,length,strand):
            strand_list=[]
            if strand=="fw":
                strand_list=[1,0]
            else:
                strand_list=[0,1]

            primer=primer3.bindings.design_primers(
                        seq_args={'SEQUENCE_TEMPLATE': sequence},
                        global_args={
                            'PRIMER_PICK_LEFT_PRIMER': "{0}".format(strand_list[0]),
                            'PRIMER_PICK_RIGHT_PRIMER': "{0}".format(strand_list[1]),
                            'PRIMER_PICK_INTERNAL_OLIGO': 0,
                            'PRIMER_OPT_SIZE': length,
                            'PRIMER_MIN_SIZE': length,
                            'PRIMER_MAX_SIZE': length,
                            'PRIMER_OPT_TM': self.opt_tm,
                            'PRIMER_MIN_TM': self.min_tm,
                            'PRIMER_MAX_TM': self.max_tm,
                            'PRIMER_OPT_GC': self.opt_GC,
                            'PRIMER_MIN_GC': self.min_GC,
                            'PRIMER_MAX_GC': self.max_GC,
                            'PRIMER_MAX_POLY_X': int(self.max_poly_x),
                            'PRIMER_MAX_NS_ACCEPTED': 0,
                            'PRIMER_MAX_SELF_ANY': int(self.max_self_any),
                            'PRIMER_MAX_SELF_END': int(self.max_self_end)
                        })
            
            primer_seq="-"
            primer_tm="-"
            primer_gc="-"
            primer_ln="-"
            if strand=="fw":
                if primer['PRIMER_LEFT_NUM_RETURNED']>0:
                    primer_seq=primer['PRIMER_LEFT_0_SEQUENCE']
                    primer_tm=int(primer['PRIMER_LEFT_0_TM']*10)/10
                    primer_gc=int(primer['PRIMER_LEFT_0_GC_PERCENT']*10)/10
                    primer_ln=primer['PRIMER_LEFT_0'][1]
            else:
                if primer['PRIMER_RIGHT_NUM_RETURNED']>0:
                    primer_seq=primer['PRIMER_RIGHT_0_SEQUENCE']
                    primer_tm=int(primer['PRIMER_RIGHT_0_TM']*10)/10
                    primer_gc=int(primer['PRIMER_RIGHT_0_GC_PERCENT']*10)/10
                    primer_ln=primer['PRIMER_RIGHT_0'][1]
            return primer_seq,primer_tm,primer_gc,primer_ln

        with open('{0}/target_SNP_selection/SNP.txt'.format(self.output_dir)) as f:
            for line in f:
                line=line.strip()
                array=line.split("\t")
                position=int(array[1])-self.s_posi

                fw_primer_list=[]
                rv_primer_list=[]

                select_s_posi=int(position-self.max_ln+1)
                select_e_posi=int(position+self.max_ln)
                if select_s_posi<0 or self.e_posi-self.s_posi+1<select_e_posi:
                    pass
                else:
                    for sample_count,sample in enumerate([self.aname,self.bname]):
                        fw_primer="-"
                        fw_tm="-"
                        fw_gc="-"
                        fw_ln="-"

                        rv_primer="-"
                        rv_tm="-"
                        rv_gc="-"
                        rv_ln="-"

                        seq="" # raw_sample_seq
                        with open('{0}/target_SNP_selection/{1}.{2}.{3}-{4}.seq'.format(self.output_dir,sample,self.chr,self.s_posi,self.e_posi)) as seq_f:
                            for whole_seq in seq_f:
                                whole_seq=whole_seq.strip()
                                seq=whole_seq[select_s_posi:select_e_posi]
                        
                        seq2="" # replace N in error
                        for count,i in enumerate(range(select_s_posi+self.s_posi,select_e_posi+self.s_posi)):
                            base=seq[count]
                            i=str(i)
                            if i in self.dict_depth_error:
                                if self.dict_depth_error[i].split(":")[sample_count]!=".":
                                    base="N"
                            if i in self.dict_variant_error:
                                variant_type=self.dict_variant_error[i].split(":")[sample_count]
                                
                                if variant_type!="." and variant_type!="s" and int(i)-self.s_posi!=position:
                                    base="N"
                            seq2=seq2+base

                        # fw作成
                        if "N" in seq2[self.max_ln-self.min_ln:self.max_ln-1]:
                            pass
                        else:
                            for i in self.len_list:
                                if fw_primer=="-":
                                    seq3=seq2[self.max_ln-i:self.max_ln] # 長さを調節
                                    base_list=[]
                                    if seq3[-3]=="A":
                                        base_list=["C","T","G"]
                                    elif seq3[-3]=="T":
                                        base_list=["C","A","G"]
                                    elif seq3[-3]=="C":
                                        base_list=["A","T","G"]
                                    elif seq3[-3]=="G":
                                        base_list=["A","T","C"]
                                    for pattern in base_list:
                                        if fw_primer=="-":
                                            seq4=seq3[0:i-3]+pattern+seq3[i-2:i] # 3'末端の2bp上流に置換
                                            if seq4[-1]==seq4[-2] and seq4[-1]==seq4[-3]:
                                                pass
                                            else:
                                                fw_primer,fw_tm,fw_gc,fw_ln=run_primer3(seq4,i,"fw")
                        fw_primer_list.append(fw_primer)
                        fw_primer_list.append(fw_tm)
                        fw_primer_list.append(fw_gc)
                        fw_primer_list.append(fw_ln)

                        # rv作成
                        if "N" in seq2[self.max_ln-1:self.max_ln+self.min_ln-1]:
                            pass
                        else:
                            for i in self.len_list:
                                if rv_primer=="-":
                                    seq3=seq2[self.max_ln-1:self.max_ln+i] # 長さを調節
                                    base_list=[]
                                    if seq3[2]=="A":
                                        base_list=["C","T","G"]
                                    elif seq3[2]=="T":
                                        base_list=["C","A","G"]
                                    elif seq3[2]=="C":
                                        base_list=["A","T","G"]
                                    elif seq3[2]=="G":
                                        base_list=["A","T","C"]
                                    for pattern in base_list:
                                        if rv_primer=="-":
                                            seq4=seq3[0:2]+pattern+seq3[3:i] # 3'末端の2bp上流に置換
                                            if seq4[0]==seq4[1] and seq4[0]==seq4[2]:
                                                pass
                                            else:
                                                rv_primer,rv_tm,rv_gc,rv_ln=run_primer3(seq4,i,"rv")
                        rv_primer_list.append(rv_primer)
                        rv_primer_list.append(rv_tm)
                        rv_primer_list.append(rv_gc)
                        rv_primer_list.append(rv_ln)

                if fw_primer_list.count('-')!=8 or rv_primer_list.count('-')!=8:
                    output_test_primers.write("{0}\t{1}\t{2}_Fw\t{3}\t{4}\t{5}\t{6}\t{7}_Fw\t{8}\t{9}\t{10}\t{11}\t{2}_Rv\t{12}\t{13}\t{14}\t{15}\t{7}_Rv\t{16}\t{17}\t{18}\t{19}\n".format(
                                                                                                            array[0],array[1],
                                                                                                            self.aname,fw_primer_list[0],fw_primer_list[1],fw_primer_list[2],fw_primer_list[3],
                                                                                                            self.bname,fw_primer_list[4],fw_primer_list[5],fw_primer_list[6],fw_primer_list[7],
                                                                                                            rv_primer_list[0],rv_primer_list[1],rv_primer_list[2],rv_primer_list[3],
                                                                                                            rv_primer_list[4],rv_primer_list[5],rv_primer_list[6],rv_primer_list[7]
                                                                                                            ))
        output_test_primers.close()   
    def run(self):
        print('[MarkMaker:{}] Start ARMS_preparation'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.preparation()
        print('[MarkMaker:{}] Make primer'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.make_primer()

        output = open('{0}/log/ARMS_preparation.log'.format(self.output_dir), "w")
        count1,count2,count3,count4 = 0,0,0,0
        with open('{0}/ARMS_preparation/made_primers.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                line=line.strip()
                array=line.split("\t")
                if array[3]!="-":
                    count1 += 1
                if array[8]!="-":
                    count2 += 1
                if array[13]!="-":
                    count3 += 1
                if array[18]!="-":
                    count4 += 1
        output.write("{0} {1} specific Fw primers were made\n".format(count1,self.aname))
        output.write("{0} {1} specific Fw primers were made\n".format(count2,self.bname))
        output.write("{0} {1} specific Rv primers were made\n".format(count3,self.aname))
        output.write("{0} {1} specific Rv primers were made\n".format(count4,self.bname))
        
        print('[MarkMaker:{}] Finish ARMS_preparation'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))