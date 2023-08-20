import os
import re
import primer3
import math
from datetime import datetime

class CAPS(object):
    def __init__(self,output_dir,restriction_enzyme,recipe,PCR_max_size,PCR_min_size,fragment_min_size,thread,hope_html):

        self.output_dir,self.restriction_enzyme,self.recipe=output_dir,restriction_enzyme,recipe
        self.PCR_max_size,self.PCR_min_size,self.fragment_min_size,self.thread=PCR_max_size,PCR_min_size,fragment_min_size,thread
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

        os.makedirs('{0}/CAPS'.format(self.output_dir), exist_ok=True)
        
        # 制限酵素情報
        self.enzyme_dict={}
        with open(self.restriction_enzyme) as f:
            for line in f:
                
                line=line.strip()
                array=line.split()
                if len(array)>=2:
                    name=array[0]
                    bases=array[1]
                    if name in self.enzyme_dict.keys():
                        self.enzyme_dict[name]="{0},{1}".format(self.enzyme_dict[name],bases)
                    else:
                        self.enzyme_dict[name]=bases

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

    def search_target_position(self):
        snp_posi_dict={}
        with open('{0}/target_SNP_selection/SNP.txt'.format(self.output_dir)) as f:
            for line in f:
                line=line.strip()
                array=line.split("\t")
                posi=int(array[1])
                snp_posi_dict[posi]=1

        for sample_name in [self.aname,self.bname]:
            output = open('{0}/CAPS/CAPS_site_cut_{1}.txt'.format(self.output_dir,sample_name), "w")
            with open('{0}/target_SNP_selection/{1}.{2}.{3}-{4}.seq'.format(self.output_dir,sample_name,self.chr,self.s_posi,self.e_posi)) as f:
                for line in f:
                    line=line.strip()
                    for enzyme_name in self.enzyme_dict.keys():
                        enzyme_base=[]
                        if "," in self.enzyme_dict[enzyme_name]:
                            enzyme_base_list=self.enzyme_dict[enzyme_name].split(",")
                            for i in enzyme_base_list:
                                enzyme_base.append(i)
                        else:
                            enzyme_base=[self.enzyme_dict[enzyme_name]]

                        for target_base in enzyme_base:
                            matches = re.finditer(target_base, line)
                            for match_site in matches:
                                clear="yes"
                                out=""
                                for i in range(match_site.start()+self.s_posi,match_site.end()+self.s_posi):
                                    if str(i) in self.dict_variant_error:
                                        if self.dict_variant_error[str(i)].split(":")[0]!="s":
                                            clear="no"
                                    if str(i) in self.dict_depth_error:
                                        clear="no"
                                    if i in snp_posi_dict.keys():
                                        out="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(self.chr,i,line[i-self.s_posi],enzyme_name,target_base,match_site.start()+self.s_posi,match_site.end()+self.s_posi-1)
                                if clear=="yes" and out != "":
                                    output.write("{0}\n".format(out))
            output.close()
        
    def make_primer(self):
        output_primers = open('{0}/CAPS/made_primers.txt'.format(self.output_dir), "w")
        for sample_name in [self.aname,self.bname]:
            with open('{0}/CAPS/CAPS_site_cut_{1}.txt'.format(self.output_dir,sample_name)) as f:
                for line in f:
                    line=line.strip()
                    array=line.split("\t")
                    select_s_posi=int(array[5])-self.s_posi-(self.PCR_max_size-self.fragment_min_size)
                    select_e_posi=int(array[6])-self.s_posi+1+(self.PCR_max_size-self.fragment_min_size)
                    seq_target_s_posi=self.PCR_max_size-self.fragment_min_size-self.fragment_min_size-len(array[5])
                    seq_target_length=len(array[4])+self.fragment_min_size+self.fragment_min_size
                    if select_s_posi<0 or self.e_posi-self.s_posi+1<select_e_posi:
                        pass
                    else:
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
                                'PRIMER_PRODUCT_SIZE_RANGE': '{0}-{1}'.format(self.PCR_min_size,self.PCR_max_size),
                                'SEQUENCE_TARGET': '{0},{1}'.format(seq_target_s_posi,seq_target_length)
                            })

                        if primer['PRIMER_PAIR_NUM_RETURNED']>0:
                            output_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n".format(
                                array[0],array[1],array[3],array[4],sample_name,select_s_posi,select_e_posi,
                                primer['PRIMER_LEFT_0_SEQUENCE'],primer['PRIMER_LEFT_0'],int(primer['PRIMER_LEFT_0_TM']*10)/10,int(primer['PRIMER_LEFT_0_GC_PERCENT']*10)/10,
                                primer['PRIMER_RIGHT_0_SEQUENCE'],primer['PRIMER_RIGHT_0'],int(primer['PRIMER_RIGHT_0_TM']*10)/10,int(primer['PRIMER_RIGHT_0_GC_PERCENT']*10)/10
                            ))
        output_primers.close()
    def make_html(self):
        with open('{0}/CAPS/made_primers.txt'.format(self.output_dir)) as marker_txt:
            for marker_line in marker_txt:
                marker_line=marker_line.strip()
                marker_array=marker_line.split("\t")
                info_chr=marker_array[0]
                info_position=int(marker_array[1])
                info_restriction_name=marker_array[2]
                info_restriction_seq=marker_array[3]
                info_Fw_seq=marker_array[7]
                info_Rv_seq=marker_array[11]
                marker_f_posi=marker_array[8]
                marker_f_posi=marker_f_posi.replace("[","")
                marker_f_posi=marker_f_posi.replace("]","")
                marker_f_posi=marker_f_posi.replace(",","")
                marker_f_posi_array=marker_f_posi.split()
                info_Fw_posi=int(marker_f_posi_array[0])
                info_Fw_len=int(marker_f_posi_array[1])
                marker_r_posi=marker_array[12]
                marker_r_posi=marker_r_posi.replace("[","")
                marker_r_posi=marker_r_posi.replace("]","")
                marker_r_posi=marker_r_posi.replace(",","")
                marker_r_posi_array=marker_r_posi.split()
                info_Rv_posi=int(marker_r_posi_array[0])
                info_Rv_len=int(marker_r_posi_array[1])
                info_Fw_tm=float(marker_array[9])
                info_Rv_tm=float(marker_array[13])
                info_Fw_gc=float(marker_array[10])
                info_Rv_gc=float(marker_array[14])
                info_seq_s_posi=int(marker_array[5])+info_Fw_posi
                info_seq_e_posi=info_seq_s_posi+info_Rv_posi-info_Fw_posi+1
                cut_sample=marker_array[4]

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

                # 制限酵素サイトをサーチ
                enzyme_base=[]
                if "," in self.enzyme_dict[info_restriction_name]:
                    enzyme_base_list=self.enzyme_dict[info_restriction_name].split(",")
                    for i in enzyme_base_list:
                        enzyme_base.append(i)
                else:
                    enzyme_base=[self.enzyme_dict[info_restriction_name]]

                a_rest_posi=[]
                for target_base in enzyme_base:
                    matches = re.finditer(target_base, aseq)
                    for match_site in matches:
                        for i in range(match_site.start(),match_site.end()):
                            a_rest_posi.append(i)
                b_rest_posi=[]
                for target_base in enzyme_base:
                    matches = re.finditer(target_base, bseq)
                    for match_site in matches:
                        for i in range(match_site.start(),match_site.end()):
                            b_rest_posi.append(i)

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

                    if count in a_rest_posi:
                        dict_a_seq[count]='<font color="#0000ff">{0}</font>'.format(dict_a_seq[count])
                    if count in b_rest_posi:
                        dict_b_seq[count]='<font color="#0000ff">{0}</font>'.format(dict_b_seq[count])

                    if count<=info_Fw_len-1:
                        avariant=">"
                        bvariant=">"
                    if count>=info_Rv_posi-info_Fw_posi+1-info_Rv_len:
                        avariant="<"
                        bvariant="<"

                    dict_a_seq2[count]=avariant
                    dict_b_seq2[count]=bvariant

                output_html = open('{0}/CAPS/html/{1}.{2}.{3}_cut.{4}.html'.format(self.output_dir,info_chr,info_position,cut_sample,info_restriction_name), "w")
                output_html.write('<!DOCTYPE html><html><head><meta charset="UTF-8">')
                output_html.write('<title>{0}.{1}.{2}_cut.{3}</title>\n'.format(info_chr,info_position,cut_sample,info_restriction_name))
                output_html.write('<style>body {font-family: Monaco, monospace;}</style></head><body><h1>CAPS marker</h1>\n')
                output_html.write('<table><tr><th>Chrmosome</th><th>Position</th></tr><tr><td>{0}</td><td>{1}-{2}</td></tr></table>\n'.format(info_chr,real_s_posi,real_e_posi))
                output_html.write('<table><tr><th>Target SNP</th></tr><tr><td>{0}</td><td>({1})</td></tr></table>\n'.format(info_position,info_position-self.s_posi-info_seq_s_posi+1))
                output_html.write('<table><tr><th>Restriction</th><th>Sequence</th><th>Sample</th></tr><tr><td>{0}</td><td>{1}</td><td>{2}</td></tr></table>\n'.format(info_restriction_name,info_restriction_seq,cut_sample))
                output_html.write('<table><tr><th>Primers</th><th>Sequence</th><th>Length</th><th>Tm</th><th>GC%</th></tr>')
                output_html.write('<tr><td>Fw</td><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>'.format(info_Fw_seq,info_Fw_len,info_Fw_tm,info_Fw_gc))
                output_html.write('<tr><td>Rv</td><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr></table>\n'.format(info_Rv_seq,info_Rv_len,info_Rv_tm,info_Rv_gc))
                output_html.write('<table><tr><th>PCR product size</th></tr><tr><td>Fw & Rv : {0}</td></tr></table>\n'.format(info_Rv_posi-info_Fw_posi+1))
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
        print('[MarkMaker:{}] Start CAPS'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.preparation()
        print('[MarkMaker:{}] Call CAPS position'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.search_target_position()
        print('[MarkMaker:{}] Make primer'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.make_primer()
        if self.hope_html=="yes":
            os.makedirs('{0}/CAPS/html'.format(self.output_dir), exist_ok=True)
            print('[MarkMaker:{}] Make HTML'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            self.make_html()

        output = open('{0}/log/CAPS.log'.format(self.output_dir), "w")
        count,count1,count2 = 0,0,0
        with open('{0}/CAPS/made_primers.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                count += 1
        with open('{0}/CAPS/CAPS_site_cut_{1}.txt'.format(self.output_dir,self.aname), 'r') as f:
            for line in f:
                count1 += 1
        with open('{0}/CAPS/CAPS_site_cut_{1}.txt'.format(self.output_dir,self.bname), 'r') as f:
            for line in f:
                count2 += 1
        output.write("{0} CAPS sites were called in {1}\n".format(count1,self.aname))
        output.write("{0} CAPS sites were called in {1}\n".format(count2,self.bname))
        output.write("{0} CAPS maekers were made\n".format(count))
        print('[MarkMaker:{}] Finish CAPS'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        