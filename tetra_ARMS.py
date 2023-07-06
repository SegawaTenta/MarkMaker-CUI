import os
import primer3
import math
from datetime import datetime

class tetra_ARMS(object):
    def __init__(self,output_dir,recipe,thread,fisrt_size_min,fisrt_size_max,second_size_min,second_size_max,hope_html):

        self.output_dir,self.recipe,self.thread=output_dir,recipe,thread
        self.fisrt_size_min,self.fisrt_size_max,self.second_size_min,self.second_size_max=fisrt_size_min,fisrt_size_max,second_size_min,second_size_max
        self.hope_html=hope_html

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
        self.max_compl_any=6
        self.max_compl_end=3

        with open(self.recipe) as f:
            for line in f:
                line=line.strip()
                array=line.split("=")
                set_name=array[0]
                info=float(array[1])
                if "PRIMER_OPT_SIZE"==set_name:
                    self.opt_ln=info
                elif "PRIMER_MIN_SIZE"==set_name:
                    self.min_ln=info
                elif "PRIMER_MAX_SIZE"==set_name:
                    self.max_ln=info
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
                elif "PRIMER_PAIR_MAX_COMPL_ANY"==set_name:
                    self.max_compl_any=info
                elif "PRIMER_PAIR_MAX_COMPL_END"==set_name:
                    self.max_compl_end=info

        os.makedirs('{0}/tetra_ARMS'.format(self.output_dir), exist_ok=True)

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

    def make_primer(self):
        output_primers = open('{0}/tetra_ARMS/made_primers.txt'.format(self.output_dir), "w")
        output_site = open('{0}/tetra_ARMS/tetra_site.txt'.format(self.output_dir), "w")
        with open('{0}/ARMS_preparation/made_primers.txt'.format(self.output_dir)) as f:
            for line in f:
                line=line.strip()
                array=line.split("\t")
                position=int(array[1])
                a_fw=array[3]
                a_rv=array[13]
                b_fw=array[8]
                b_rv=array[18]

                select_s_posi=position-self.s_posi-self.fisrt_size_max
                select_e_posi=position-self.s_posi+self.second_size_max
                seq_target_s_posi=self.fisrt_size_max-self.fisrt_size_min
                seq_target_length=self.fisrt_size_min+1+self.second_size_min
                
                can_make="no"
                if select_s_posi<0 or self.e_posi-self.s_posi+1<select_e_posi:
                        pass
                else:
                    if a_fw!="-" and b_rv!="-":
                        can_make="yes"
                        a_fw_name=array[2]
                        a_fw_tm=float(array[4])
                        a_fw_gc=float(array[5])
                        a_fw_ln=int(array[6])
                        b_rv_name=array[17]
                        b_rv_tm=float(array[19])
                        b_rv_gc=float(array[20])
                        b_rv_ln=int(array[21])
                        seq=""
                        with open('{0}/target_SNP_selection/shared.{1}.{2}-{3}.seq'.format(self.output_dir,self.chr,self.s_posi,self.e_posi)) as seq_f:
                            for whole_seq in seq_f:
                                whole_seq=whole_seq.strip()
                                seq=whole_seq[select_s_posi:select_e_posi]
                        
                        primer=primer3.bindings.design_primers(
                            seq_args={'SEQUENCE_TEMPLATE': seq},
                            global_args={
                                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                                'PRIMER_OPT_SIZE': int(self.opt_ln),
                                'PRIMER_MIN_SIZE': int(self.min_ln),
                                'PRIMER_MAX_SIZE': int(self.max_ln),
                                'PRIMER_OPT_TM': self.opt_tm,
                                'PRIMER_MIN_TM': self.min_tm,
                                'PRIMER_MAX_TM': self.max_tm,
                                'PRIMER_OPT_GC': self.opt_GC,
                                'PRIMER_MIN_GC': self.min_GC,
                                'PRIMER_MAX_GC': self.max_GC,
                                'PRIMER_MAX_POLY_X': int(self.max_poly_x),
                                'PRIMER_MAX_NS_ACCEPTED': 0,
                                'PRIMER_MAX_SELF_ANY': int(self.max_self_any),
                                'PRIMER_MAX_SELF_END': int(self.max_self_end),
                                'PRIMER_PAIR_MAX_COMPL_ANY': int(self.max_compl_any),
                                'PRIMER_PAIR_MAX_COMPL_END': int(self.max_compl_end),
                                'PRIMER_PRODUCT_SIZE_RANGE': '{0}-{1}'.format(self.fisrt_size_min+self.second_size_min+1,self.fisrt_size_max+self.second_size_max+1),
                                'SEQUENCE_TARGET': '{0},{1}'.format(seq_target_s_posi,seq_target_length)
                            })
                        
                        # print(seq,select_s_posi,select_e_posi)

                        if primer['PRIMER_PAIR_NUM_RETURNED']>0:
                            output_primers.write("{0}\t{1}\t{2}\t{3}\t".format(array[0],position,select_s_posi,select_e_posi))
                            output_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(a_fw_name,a_fw,a_fw_ln,a_fw_tm,a_fw_gc))
                            output_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(b_rv_name,b_rv,b_rv_ln,b_rv_tm,b_rv_gc))
                            output_primers.write("Fw\t{0}\t{1}\t{2}\t{3}\t".format(primer['PRIMER_LEFT_0_SEQUENCE'],primer['PRIMER_LEFT_0'],int(primer['PRIMER_LEFT_0_TM']*10)/10,int(primer['PRIMER_LEFT_0_GC_PERCENT']*10)/10))
                            output_primers.write("Rv\t{0}\t{1}\t{2}\t{3}\n".format(primer['PRIMER_RIGHT_0_SEQUENCE'],primer['PRIMER_RIGHT_0'],int(primer['PRIMER_RIGHT_0_TM']*10)/10,int(primer['PRIMER_RIGHT_0_GC_PERCENT']*10)/10))
                            can_make="did"

                    if a_rv!="-" and b_fw!="-" and can_make!="did":
                        can_make="yes"
                        b_fw_name=array[7]
                        b_fw_tm=float(array[9])
                        b_fw_gc=float(array[10])
                        b_fw_ln=int(array[11])
                        a_rv_name=array[12]
                        a_rv_tm=float(array[14])
                        a_rv_gc=float(array[15])
                        a_rv_ln=int(array[16])
                        seq=""
                        with open('{0}/target_SNP_selection/shared.{1}.{2}-{3}.seq'.format(self.output_dir,self.chr,self.s_posi,self.e_posi)) as seq_f:
                            for whole_seq in seq_f:
                                whole_seq=whole_seq.strip()
                                seq=whole_seq[select_s_posi:select_e_posi]
                        
                        primer=primer3.bindings.design_primers(
                            seq_args={'SEQUENCE_TEMPLATE': seq},
                            global_args={
                                'PRIMER_PICK_INTERNAL_OLIGO': 0,
                                'PRIMER_OPT_SIZE': int(self.opt_ln),
                                'PRIMER_MIN_SIZE': int(self.min_ln),
                                'PRIMER_MAX_SIZE': int(self.max_ln),
                                'PRIMER_OPT_TM': self.opt_tm,
                                'PRIMER_MIN_TM': self.min_tm,
                                'PRIMER_MAX_TM': self.max_tm,
                                'PRIMER_OPT_GC': self.opt_GC,
                                'PRIMER_MIN_GC': self.min_GC,
                                'PRIMER_MAX_GC': self.max_GC,
                                'PRIMER_MAX_POLY_X': int(self.max_poly_x),
                                'PRIMER_MAX_NS_ACCEPTED': 0,
                                'PRIMER_MAX_SELF_ANY': int(self.max_self_any),
                                'PRIMER_MAX_SELF_END': int(self.max_self_end),
                                'PRIMER_PAIR_MAX_COMPL_ANY': int(self.max_compl_any),
                                'PRIMER_PAIR_MAX_COMPL_END': int(self.max_compl_end),
                                'PRIMER_PRODUCT_SIZE_RANGE': '{0}-{1}'.format(self.fisrt_size_min+self.second_size_min+1,self.fisrt_size_max+self.second_size_max+1),
                                'SEQUENCE_TARGET': '{0},{1}'.format(seq_target_s_posi,seq_target_length)
                            })

                        if primer['PRIMER_PAIR_NUM_RETURNED']>0:
                            output_primers.write("{0}\t{1}\t{2}\t{3}\t".format(array[0],position,select_s_posi,select_e_posi))
                            output_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(b_fw_name,b_fw,b_fw_ln,b_fw_tm,b_fw_gc))
                            output_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(a_rv_name,a_rv,a_rv_ln,a_rv_tm,a_rv_gc))
                            output_primers.write("Fw\t{0}\t{1}\t{2}\t{3}\t".format(primer['PRIMER_LEFT_0_SEQUENCE'],primer['PRIMER_LEFT_0'],int(primer['PRIMER_LEFT_0_TM']*10)/10,int(primer['PRIMER_LEFT_0_GC_PERCENT']*10)/10))
                            output_primers.write("Rv\t{0}\t{1}\t{2}\t{3}\n".format(primer['PRIMER_RIGHT_0_SEQUENCE'],primer['PRIMER_RIGHT_0'],int(primer['PRIMER_RIGHT_0_TM']*10)/10,int(primer['PRIMER_RIGHT_0_GC_PERCENT']*10)/10))

                if can_make=="yes" or can_make=="did":
                    output_site.write("{0}\n".format(line))
    
        output_primers.close()
        output_site.close()
    
    def make_html(self):
        with open('{0}/tetra_ARMS/made_primers.txt'.format(self.output_dir)) as marker_txt:
            for marker_line in marker_txt:
                marker_line=marker_line.strip()
                marker_array=marker_line.split("\t")
                info_chr=marker_array[0]
                info_position=int(marker_array[1])
                info_Fw_seq=marker_array[15]
                info_Rv_seq=marker_array[20]
                info_s_Fw_seq=marker_array[5]
                info_s_Rv_seq=marker_array[10]

                marker_f_posi=marker_array[16]
                marker_f_posi=marker_f_posi.replace("[","")
                marker_f_posi=marker_f_posi.replace("]","")
                marker_f_posi=marker_f_posi.replace(",","")
                marker_f_posi_array=marker_f_posi.split()
                info_Fw_posi=int(marker_f_posi_array[0])
                info_Fw_len=int(marker_f_posi_array[1])

                marker_r_posi=marker_array[21]
                marker_r_posi=marker_r_posi.replace("[","")
                marker_r_posi=marker_r_posi.replace("]","")
                marker_r_posi=marker_r_posi.replace(",","")
                marker_r_posi_array=marker_r_posi.split()
                info_Rv_posi=int(marker_r_posi_array[0])
                info_Rv_len=int(marker_r_posi_array[1])

                info_Fw_tm=float(marker_array[17])
                info_Rv_tm=float(marker_array[22])
                info_s_Fw_tm=float(marker_array[7])
                info_s_Rv_tm=float(marker_array[12])

                info_Fw_gc=float(marker_array[18])
                info_Rv_gc=float(marker_array[23])
                info_s_Fw_gc=float(marker_array[8])
                info_s_Rv_gc=float(marker_array[13])

                info_seq_s_posi=int(marker_array[2])+info_Fw_posi
                info_seq_e_posi=info_seq_s_posi+info_Rv_posi-info_Fw_posi+1

                # 配列取得
                aseq=""
                bseq=""
                with open('{0}/target_SNP_selection/{1}.{2}.{3}-{4}.seq'.format(self.output_dir,self.aname,self.chr,self.s_posi,self.e_posi)) as seq_f:
                    for whole_seq in seq_f:
                        whole_seq=whole_seq.strip()
                        aseq=whole_seq[info_seq_s_posi:info_seq_e_posi]
                with open('{0}/target_SNP_selection/{1}.{2}.{3}-{4}.seq'.format(self.output_dir,self.bname,self.chr,self.s_posi,self.e_posi)) as seq_f:
                    for whole_seq in seq_f:
                        whole_seq=whole_seq.strip()
                        bseq=whole_seq[info_seq_s_posi:info_seq_e_posi]
                dict_a_seq={}
                dict_b_seq={}
                for i in range(len(aseq)):
                    dict_a_seq[i]=aseq[i]
                for i in range(len(bseq)):
                    dict_b_seq[i]=bseq[i]

                real_s_posi=info_seq_s_posi+self.s_posi
                real_e_posi=info_seq_e_posi+self.s_posi

                dict_a_seq2={}
                dict_b_seq2={}
                for count,i in enumerate(range(real_s_posi,real_e_posi)):
                    # avariant="&nbsp;"
                    avariant="&nbsp;"
                    bvariant="&nbsp;"

                    if str(i) in self.dict_depth_error.keys():
                        variant_array=self.dict_depth_error[str(i)].split(":")
                        a_depth=variant_array[0]
                        b_depth=variant_array[1]
                        
                        if a_depth=="0":
                                avariant="0"
                        if a_depth=="L":
                                avariant="L"
                        if a_depth=="H":
                                avariant="H"
                        
                        if b_depth=="0":
                                bvariant="0"
                        if b_depth=="L":
                                bvariant="L"
                        if b_depth=="H":
                                bvariant="H"

                    if str(i) in self.dict_variant_error.keys():
                        variant_array=self.dict_variant_error[str(i)].split(":")
                        a_variant=variant_array[0]
                        b_variant=variant_array[1]
                        if a_variant=="+" or a_variant=="-" or a_variant=="M":
                            avariant=a_variant
                            dict_a_seq[count]='<font color="#008000">{0}</font>'.format(dict_a_seq[count])
                        elif a_variant=="m" or a_variant=="v":
                            dict_a_seq[count]='<font color="#008000">{0}</font>'.format(dict_a_seq[count])
                        elif a_variant=="s":
                            dict_a_seq[count]='<font color="#ff0000">{0}</font>'.format(dict_a_seq[count])

                        if b_variant=="+" or b_variant=="-" or b_variant=="M":
                            bvariant=b_variant
                            dict_b_seq[count]='<font color="#008000">{0}</font>'.format(dict_b_seq[count])
                        elif b_variant=="m" or b_variant=="v":
                            dict_b_seq[count]='<font color="#008000">{0}</font>'.format(dict_b_seq[count])
                        elif b_variant=="s":
                            dict_b_seq[count]='<font color="#ff0000">{0}</font>'.format(dict_b_seq[count])

                    if count<=info_Fw_len-1:
                        avariant=">"
                        bvariant=">"
                    if count>=info_Rv_posi-info_Fw_posi+1-info_Rv_len:
                        avariant="<"
                        bvariant="<"
                    if i==info_position:
                        avariant="*"
                        bvariant="*"
                    elif i>=info_position-len(info_s_Fw_seq)+1 and i<info_position:
                        avariant='<font color="#a260bf">></font>'
                        bvariant='<font color="#a260bf">></font>'
                    elif i>info_position and i<=info_position+len(info_s_Rv_seq)-1:
                        avariant='<font color="#fd7e00"><</font>'
                        bvariant='<font color="#fd7e00"><</font>'

                    dict_a_seq2[count]=avariant
                    dict_b_seq2[count]=bvariant

                output_html = open('{0}/tetra_ARMS/html/{1}.{2}.html'.format(self.output_dir,info_chr,info_position), "w")
                output_html.write('<!DOCTYPE html><html><head><meta charset="UTF-8">')
                output_html.write('<title>{0}.{1}</title>\n'.format(info_chr,info_position))
                output_html.write('<style>body {font-family: Monaco, monospace;}</style></head><body><h1>tetra ARMS marker</h1>\n')
                output_html.write('<table><tr><th>Chrmosome</th><th>Position</th></tr><tr><td>{0}</td><td>{1}-{2}</td></tr></table>\n'.format(info_chr,real_s_posi,real_e_posi))
                output_html.write('<table><tr><th>Target SNP</th></tr><tr><td>{0}</td><td>({1})</td></tr></table>\n'.format(info_position,info_position-self.s_posi-info_seq_s_posi+1))
                output_html.write('<table><tr><th>Primers</th><th>Sequence</th><th>Length</th><th>Tm</th><th>GC%</th></tr>')
                output_html.write('<tr><td>Fw</td><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>'.format(info_Fw_seq,info_Fw_len,info_Fw_tm,info_Fw_gc))
                output_html.write('<tr><td>Rv</td><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>'.format(info_Rv_seq,info_Rv_len,info_Rv_tm,info_Rv_gc))
                output_html.write('<tr><td>{5}{0}{6}</td><td>{5}{1}{6}</td><td>{5}{2}{6}</td><td>{5}{3}{6}</td><td>{5}{4}{6}</td></tr>'.format(marker_array[4],info_s_Fw_seq,len(info_s_Fw_seq),info_s_Fw_tm,info_s_Fw_gc,'<font color="#a260bf">','</font>'))
                output_html.write('<tr><td>{5}{0}{6}</td><td>{5}{1}{6}</td><td>{5}{2}{6}</td><td>{5}{3}{6}</td><td>{5}{4}{6}</td></tr></font></table>\n'.format(marker_array[9],info_s_Rv_seq,len(info_s_Rv_seq),info_s_Rv_tm,info_s_Rv_gc,'<font color="#fd7e00">','</font>'))
                output_html.write('<table><tr><th>PCR product size</th></tr><tr><td>{0}</td></tr></table>\n'.format(info_Rv_posi-info_Fw_posi+1))
                output_html.write('<h2>{0}</h2>\n'.format(self.aname))

                divid_num=math.ceil(len(aseq)/100)
                for i in range(divid_num):
                    s_num=i*100
                    e_num=s_num+100
                    output_html.write('<p style="line-height: 0.2";>')
                    for count in range(s_num,e_num):
                        if count in dict_a_seq:
                            output_html.write('{0}'.format(dict_a_seq[count]))
                    output_html.write('</p>\n')
                    output_html.write('<p style="line-height: 0.2";>')
                    for count in range(s_num,e_num):
                        if count in dict_a_seq2:
                            output_html.write('{0}'.format(dict_a_seq2[count]))
                    output_html.write('</p>\n')

                output_html.write('<h2>{0}</h2>\n'.format(self.bname))
                for i in range(divid_num):
                    s_num=i*100
                    e_num=s_num+100
                    output_html.write('<p style="line-height: 0.2";>')
                    for count in range(s_num,e_num):
                        if count in dict_b_seq:
                            output_html.write('{0}'.format(dict_b_seq[count]))
                    output_html.write('</p>\n')
                    output_html.write('<p style="line-height: 0.2";>')
                    for count in range(s_num,e_num):
                        if count in dict_b_seq2:
                            output_html.write('{0}'.format(dict_b_seq2[count]))
                    output_html.write('</p>\n')

                output_html.close()

    def run(self):
        print('[MarkMaker:{}] Start tetra_ARMS'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.preparation()
        print('[MarkMaker:{}] Make primer'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.make_primer()
        if self.hope_html=="yes":
            os.makedirs('{0}/tetra_ARMS/html'.format(self.output_dir), exist_ok=True)
            print('[MarkMaker:{}] Make HTML'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            self.make_html()

        output = open('{0}/log/tetra_ARMS.log'.format(self.output_dir), "w")
        count1,count2 = 0,0
        with open('{0}/tetra_ARMS/tetra_site.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                count1 += 1
        with open('{0}/tetra_ARMS/made_primers.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                count2 += 1
        output.write("{0} SNPs were called\n".format(count1))
        output.write("{0} tetra ARMS markers were made\n".format(count2))
        print('[MarkMaker:{}] Finish tetra_ARMS'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))