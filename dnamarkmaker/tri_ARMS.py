import os
import primer3
import math
from datetime import datetime

class tri_ARMS(object):
    def __init__(self,output_dir,recipe,thread,PCR_max_size,PCR_min_size,SNP_dist_min,SNP_dist_max,hope_html):

        self.output_dir,self.recipe,self.thread=output_dir,recipe,thread
        self.PCR_max_size,self.PCR_min_size,self.SNP_dist_min,self.SNP_dist_max=PCR_max_size,PCR_min_size,SNP_dist_min,SNP_dist_max
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

        os.makedirs('{0}/tri_ARMS'.format(self.output_dir), exist_ok=True)

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

    def serch_target_position(self):
        output_fw = open('{0}/tri_ARMS/tri_site_sFw.txt'.format(self.output_dir), "w")
        output_rv = open('{0}/tri_ARMS/tri_site_sRv.txt'.format(self.output_dir), "w")
        snp_posi_Afw_dict={}
        snp_posi_Bfw_dict={}
        snp_posi_Arv_dict={}
        snp_posi_Brv_dict={}
        with open('{0}/ARMS_preparation/made_primers.txt'.format(self.output_dir)) as f:
            for line in f:
                line=line.strip()
                array=line.split("\t")
                posi=int(array[1])
                a_s_fw=array[3]
                b_s_fw=array[8]
                a_s_rv=array[13]
                b_s_rv=array[18]
                
                if a_s_fw!="-":
                    for i in snp_posi_Bfw_dict.keys():
                        if posi-i>=self.SNP_dist_min and posi-i<=self.SNP_dist_max:
                            pre_array=snp_posi_Bfw_dict[i].split(":")
                            output_fw.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(array[0],i,pre_array[0],pre_array[1],pre_array[2],pre_array[3],pre_array[4],array[1],array[2],array[3],array[4],array[5],array[6]))
                    snp_posi_Afw_dict[posi]="{0}:{1}:{2}:{3}:{4}".format(array[2],array[3],array[4],array[5],array[6])

                if b_s_fw!="-":
                    for i in snp_posi_Afw_dict.keys():
                        if posi-i>=self.SNP_dist_min and posi-i<=self.SNP_dist_max:
                            pre_array=snp_posi_Afw_dict[i].split(":")
                            output_fw.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(array[0],i,pre_array[0],pre_array[1],pre_array[2],pre_array[3],pre_array[4],array[1],array[7],array[8],array[9],array[10],array[11]))
                    snp_posi_Bfw_dict[posi]="{0}:{1}:{2}:{3}:{4}".format(array[7],array[8],array[9],array[10],array[11])

                if a_s_rv!="-":
                    for i in snp_posi_Brv_dict.keys():
                        if posi-i>=self.SNP_dist_min and posi-i<=self.SNP_dist_max:
                            pre_array=snp_posi_Brv_dict[i].split(":")
                            output_rv.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(array[0],i,pre_array[0],pre_array[1],pre_array[2],pre_array[3],pre_array[4],array[1],array[12],array[13],array[14],array[15],array[16]))
                    snp_posi_Arv_dict[posi]="{0}:{1}:{2}:{3}:{4}".format(array[12],array[13],array[14],array[15],array[16])

                if b_s_rv!="-":
                    for i in snp_posi_Arv_dict.keys():
                        if posi-i>=self.SNP_dist_min and posi-i<=self.SNP_dist_max:
                            pre_array=snp_posi_Arv_dict[i].split(":")
                            output_rv.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(array[0],i,pre_array[0],pre_array[1],pre_array[2],pre_array[3],pre_array[4],array[1],array[17],array[18],array[19],array[20],array[21]))
                    snp_posi_Brv_dict[posi]="{0}:{1}:{2}:{3}:{4}".format(array[17],array[18],array[19],array[20],array[21])

    def make_primer(self):
        output_test_primers = open('{0}/tri_ARMS/made_primers.txt'.format(self.output_dir), "w")
        def run_primer3(sequence,strand):
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
                            'PRIMER_OPT_SIZE': self.opt_ln,
                            'PRIMER_MIN_SIZE': self.min_ln,
                            'PRIMER_MAX_SIZE': self.max_ln,
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
            primer_posi="-"
            if strand=="fw":
                if primer['PRIMER_LEFT_NUM_RETURNED']>0:
                    primer_seq=primer['PRIMER_LEFT_0_SEQUENCE']
                    primer_tm=int(primer['PRIMER_LEFT_0_TM']*10)/10
                    primer_gc=int(primer['PRIMER_LEFT_0_GC_PERCENT']*10)/10
                    primer_posi=primer['PRIMER_LEFT_0']
            else:
                if primer['PRIMER_RIGHT_NUM_RETURNED']>0:
                    primer_seq=primer['PRIMER_RIGHT_0_SEQUENCE']
                    primer_tm=int(primer['PRIMER_RIGHT_0_TM']*10)/10
                    primer_gc=int(primer['PRIMER_RIGHT_0_GC_PERCENT']*10)/10
                    primer_posi=primer['PRIMER_RIGHT_0']
            return primer_seq,primer_tm,primer_gc,primer_posi

        with open('{0}/tri_ARMS/tri_site_sFw.txt'.format(self.output_dir)) as f:
            for line in f:
                line=line.strip()
                array=line.split("\t")
                position=int(array[7])
                select_s_posi=position-self.s_posi+self.PCR_min_size
                select_e_posi=position-self.s_posi+self.PCR_max_size
                select_s_range=(int(array[1])-int(array[6]))-self.s_posi+1
                select_e_range=0
                if select_s_posi<0 or self.e_posi-self.s_posi+1<select_e_posi:
                        pass
                else:
                    seq=""
                    with open('{0}/target_SNP_selection/shared.{1}.{2}-{3}.seq'.format(self.output_dir,self.chr,self.s_posi,self.e_posi)) as seq_f:
                        for whole_seq in seq_f:
                            whole_seq=whole_seq.strip()
                            seq=whole_seq[select_s_posi:select_e_posi]
                    primer_seq,primer_tm,primer_gc,primer_posi=run_primer3(seq,"rv")
                    if primer_seq!="-":
                        select_e_range=select_s_posi+primer_posi[0]+1

                        product_size1=select_e_range-select_s_range
                        product_size2=select_e_range-((int(array[7])-int(array[12]))-self.s_posi+1)

                        output_test_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t".format(array[0],array[1],array[2],array[3],array[4],array[5],array[6]))
                        output_test_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(array[7],array[8],array[9],array[10],array[11],array[12]))
                        output_test_primers.write("Rv\t{0}\t{1}\t{2}\t{3}\t".format(primer_seq,primer_tm,primer_gc,primer_posi[1]))
                        output_test_primers.write("{0}\t{1}\t".format(select_s_range,select_e_range))
                        output_test_primers.write("{0}\t{1}\n".format(product_size1,product_size2))
        
        with open('{0}/tri_ARMS/tri_site_sRv.txt'.format(self.output_dir)) as f:
            for line in f:
                line=line.strip()
                array=line.split("\t")
                position=int(array[1])
                select_s_posi=position-self.s_posi-self.PCR_max_size
                select_e_posi=position-self.s_posi-self.PCR_min_size
                select_s_range=0
                select_e_range=(int(array[7])+int(array[12]))-self.s_posi
                if select_s_posi<0 or self.e_posi-self.s_posi+1<select_e_posi:
                        pass
                else:
                    seq=""
                    with open('{0}/target_SNP_selection/shared.{1}.{2}-{3}.seq'.format(self.output_dir,self.chr,self.s_posi,self.e_posi)) as seq_f:
                        for whole_seq in seq_f:
                            whole_seq=whole_seq.strip()
                            seq=whole_seq[select_s_posi:select_e_posi]
                    primer_seq,primer_tm,primer_gc,primer_posi=run_primer3(seq,"fw")
                    if primer_seq!="-":
                        select_s_range=select_s_posi+primer_posi[0]

                        product_size1=select_e_range-select_s_range
                        product_size2=((int(array[1])+int(array[6]))-self.s_posi)-select_s_range

                        output_test_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t".format(array[0],array[1],array[2],array[3],array[4],array[5],array[6]))
                        output_test_primers.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(array[7],array[8],array[9],array[10],array[11],array[12]))
                        output_test_primers.write("Fw\t{0}\t{1}\t{2}\t{3}\t".format(primer_seq,primer_tm,primer_gc,primer_posi[1]))
                        output_test_primers.write("{0}\t{1}\t".format(select_s_range,select_e_range))
                        output_test_primers.write("{0}\t{1}\n".format(product_size1,product_size2))

        output_test_primers.close()

    def make_html(self):
        def color_purple(strand):
            return '<font color="#a260bf">{0}</font>'.format(strand)
        def color_orange(strand):
            return '<font color="#fd7e00">{0}</font>'.format(strand)
        
        with open('{0}/tri_ARMS/made_primers.txt'.format(self.output_dir)) as marker_txt:
            for marker_line in marker_txt:
                marker_line=marker_line.strip()
                marker_array=marker_line.split("\t")
                info_chr=marker_array[0]
                info_share_seq=marker_array[14]
                info_s1_seq=marker_array[3]
                info_s2_seq=marker_array[9]
                info_s1_posi=int(marker_array[1])
                info_s2_posi=int(marker_array[7])
                info_share_len=int(marker_array[17])
                info_s1_len=int(marker_array[6])
                info_s2_len=int(marker_array[12])
                info_share_tm=float(marker_array[15])
                info_s1_tm=float(marker_array[4])
                info_s2_tm=float(marker_array[10])
                info_share_gc=float(marker_array[16])
                info_s1_gc=float(marker_array[5])
                info_s2_gc=float(marker_array[11])
                info_seq_s_posi=int(marker_array[18])
                info_seq_e_posi=int(marker_array[19])

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

                    if marker_array[13]=="Fw":
                        if count<=info_share_len-1:
                            avariant=">"
                            bvariant=">"
                        if i>=info_s1_posi and i<=info_s1_posi+info_s1_len-1:
                            avariant=color_purple("<")
                            bvariant=color_purple("<")
                        if i>=info_s2_posi and i<=info_s2_posi+info_s2_len-1:
                            avariant=color_orange("<")
                            bvariant=color_orange("<")

                    if marker_array[13]=="Rv":
                        if i>=real_e_posi-info_share_len:
                            avariant="<"
                            bvariant="<"
                        if i>=info_s1_posi-info_s1_len+1 and i<=info_s1_posi:
                            avariant=color_purple(">")
                            bvariant=color_purple(">")
                        if i>=info_s2_posi-info_s2_len+1 and i<=info_s2_posi:
                            avariant=color_orange(">")
                            bvariant=color_orange(">")

                    dict_a_seq2[count]=avariant
                    dict_b_seq2[count]=bvariant

                output_html=""
                if marker_array[13]=="Fw":
                    output_html = open('{0}/tri_ARMS/html/{1}.{2}.{3}.sRv.html'.format(self.output_dir,info_chr,info_s1_posi,info_s2_posi), "w")
                    output_html.write('<!DOCTYPE html><html><head><meta charset="UTF-8">')
                    output_html.write('<title>{0}.{1}.{2}.sRv</title>\n'.format(info_chr,info_s1_posi,info_s2_posi))
                else:
                    output_html = open('{0}/tri_ARMS/html/{1}.{2}.{3}.sFw.html'.format(self.output_dir,info_chr,info_s1_posi,info_s2_posi), "w")
                    output_html.write('<!DOCTYPE html><html><head><meta charset="UTF-8">')
                    output_html.write('<title>{0}.{1}.{2}.sFw</title>\n'.format(info_chr,info_s1_posi,info_s2_posi))
                output_html.write('<style>body {font-family: Monaco, monospace;}</style></head><body><h1>tri ARMS marker</h1>\n')
                output_html.write('<table><tr><th>Chrmosome</th><th>Position</th></tr><tr><td>{0}</td><td>{1}-{2}</td></tr></table>\n'.format(info_chr,real_s_posi,real_e_posi))
                output_html.write('<table><tr><th>Target SNP</th></tr><tr><td>{0}</td><td>({1})</td><td>{2}</td><td>({3})</td></tr></table>\n'.format(info_s1_posi,info_s1_posi-self.s_posi-info_seq_s_posi+1,info_s2_posi,info_s2_posi-self.s_posi-info_seq_s_posi+1))
                output_html.write('<table><tr><th>Primers</th><th>Sequence</th><th>Length</th><th>Tm</th><th>GC%</th></tr>')
                
                if marker_array[13]=="Fw":
                    output_html.write('<tr><td>Fw</td><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>'.format(info_share_seq,info_share_len,info_share_tm,info_share_gc))
                else:
                    output_html.write('<tr><td>Rv</td><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td></tr>'.format(info_share_seq,info_share_len,info_share_tm,info_share_gc))
                
                output_html.write('<tr><td>{5}{0}{6}</td><td>{5}{1}{6}</td><td>{5}{2}{6}</td><td>{5}{3}{6}</td><td>{5}{4}{6}</td></tr>'.format(marker_array[2],info_s1_seq,info_s1_len,info_s1_tm,info_s1_gc,'<font color="#a260bf">','</font>'))
                output_html.write('<tr><td>{5}{0}{6}</td><td>{5}{1}{6}</td><td>{5}{2}{6}</td><td>{5}{3}{6}</td><td>{5}{4}{6}</td></tr></font></table>\n'.format(marker_array[8],info_s2_seq,info_s2_len,info_s2_tm,info_s2_gc,'<font color="#fd7e00">','</font>'))
                
                if marker_array[13]=="Fw":
                    output_html.write('<table><tr><th>PCR product size</th></tr><tr><td>Fw & {0}{1}{2} : {3}</td></tr><tr><td>Fw & {4}{5}{2} : {6}</td></tr></table>\n'.format('<font color="#a260bf">',marker_array[2],'</font>',int(marker_array[21]),'<font color="#fd7e00">',marker_array[8],int(marker_array[20]),))
                else:
                    output_html.write('<table><tr><th>PCR product size</th></tr><tr><td>{0}{1}{2} & Rv : {3}</td></tr><tr><td> {4}{5}{2} & Rv : {6}</td></tr></table>\n'.format('<font color="#a260bf">',marker_array[2],'</font>',int(marker_array[20]),'<font color="#fd7e00">',marker_array[8],int(marker_array[21]),))
                
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
        print('[MarkMaker:{}] Start tri_ARMS'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.preparation()
        print('[MarkMaker:{}] Call target SNPs pair'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.serch_target_position()
        print('[MarkMaker:{}] Make primer'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.make_primer()
        if self.hope_html=="yes":
            print('[MarkMaker:{}] Make HTML'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            os.makedirs('{0}/tri_ARMS/html'.format(self.output_dir), exist_ok=True)
            self.make_html()
        
        output = open('{0}/log/tri_ARMS.log'.format(self.output_dir), "w")
        count1,count2,count3 = 0,0,0
        with open('{0}/tri_ARMS/tri_site_sFw.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                count1 += 1
        with open('{0}/tri_ARMS/tri_site_sRv.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                count2 += 1
        with open('{0}/tri_ARMS/made_primers.txt'.format(self.output_dir), 'r') as f:
            for line in f:
                count3 += 1
        output.write("{0} combinations of specific Fw primer were called\n".format(count1))
        output.write("{0} combinations of specific Rv primer were called\n".format(count2))
        output.write("{0} tri ARMS markers were made\n".format(count3))
        print('[MarkMaker:{}] Finish tri_ARMS'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        
