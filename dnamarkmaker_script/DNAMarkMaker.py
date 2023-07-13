#!/usr/bin/env python3

import os
import argparse
from dnamarkmaker_script.__init__ import __version__
# from datetime import datetime
# import tracemalloc

class DNAMarkMaker(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='DNAMarkMaker version {}'.format(__version__))
        self.parser.usage=('DNAMarkMaker -w <target_SNP_selection/ARMS_preparation/tri_ARMS/tetra_ARMS/CAPS>')
        self.parser.add_argument('-w','--work',
                            required=True,
                            type=str,
                            choices=['target_SNP_selection', 'CAPS', 'ARMS_preparation', 'tri_ARMS', 'tetra_ARMS'],
                            help='Choose between target_SNP_selection, CAPS, ARMS_preparation, tri_ARMS and tetra_ARMS')
        
        self.parser.add_argument('-Abam','--Abam',help='Full path of A bam')
        self.parser.add_argument('-Bbam','--Bbam',help='Full path of B bam')
        self.parser.add_argument('-Cbam','--Cbam',help='Full path of C bam')
        self.parser.add_argument('-Aname','--Aname',help='A name (A)')
        self.parser.add_argument('-Bname','--Bname',help='B name (B)')
        self.parser.add_argument('-Csim','--Csim',help='C simulation file')
        self.parser.add_argument('-reference','--reference',help='Full path of reference fasta')
        self.parser.add_argument('-position','--position',help='Target chromosome position [chr:start:end]')
        self.parser.add_argument('-o','--output_dir',help='Output directory')
        self.parser.add_argument('-min_depth','--min_depth',help='Minimum depth of target SNP (10)')
        self.parser.add_argument('-max_depth','--max_depth',help='Maximum depth of target SNP (99)')
        self.parser.add_argument('-Bhetero','--Bhetero',help='Whether to target heterozygous SNP in B (no)')
        self.parser.add_argument('-Bsim','--Bsim',help='B simulation file')
        self.parser.add_argument('-minMQ','--minMQ',help='Minimum mapping quality detected from bam (0)')
        self.parser.add_argument('-minBQ','--minBQ',help='Minimum base quality detected from bam (13)')
        self.parser.add_argument('-restriction_enzyme','--restriction_enzyme',help='Full path of restriction enzyme file')
        self.parser.add_argument('-recipe','--recipe',help='Full path of primer recipe file')
        self.parser.add_argument('-PCR_max_size','--PCR_max_size',help='Maximum size of PCR product (1000 or 700)')
        self.parser.add_argument('-PCR_min_size','--PCR_min_size',help='Minimum size of PCR product (500 or 100)')
        self.parser.add_argument('-fragment_min_size','--fragment_min_size',help='Minimum fragment size of restricted PCR product (200)')
        # self.parser.add_argument('-t','--thread',help='')
        self.parser.add_argument('-first_size','--first_size',help='Size of first band (100:500)')
        self.parser.add_argument('-second_size','--second_size',help='Size of second band (600:1000)')
        self.parser.add_argument('-SNP_dist','--SNP_dist',help='Target SNP distance (100:300)')
        self.parser.add_argument('-make_html','--make_html',help='Whether to html file (yes)')
        
        args = self.parser.parse_args()
        self.work=args.work
        self.Abam=args.Abam
        self.Bbam=args.Bbam
        self.F1bam=args.Cbam
        self.Aname=args.Aname
        self.Bname=args.Bname
        self.F1sim=args.Csim
        self.reference=args.reference
        self.position=args.position
        self.output_dir=args.output_dir
        self.min_depth=args.min_depth
        self.max_depth=args.max_depth
        self.Bhetero=args.Bhetero
        self.Bsim=args.Bsim
        self.minMQ=args.minMQ
        self.minBQ=args.minBQ
        self.restriction_enzyme=args.restriction_enzyme
        self.recipe=args.recipe
        self.PCR_max_size=args.PCR_max_size
        self.PCR_min_size=args.PCR_min_size
        self.fragment_min_size=args.fragment_min_size
        # self.thread=args.thread
        self.thread=1
        self.first_size=args.first_size
        self.second_size=args.second_size
        self.SNP_dist=args.SNP_dist
        self.make_html=args.make_html
        
        if self.work=="target_SNP_selection":
            if not args.Abam:
                self.parser.error('-Abam is required with -w target_SNP_selection')
            elif not os.path.isfile(args.Abam):
                self.parser.error('No such file {0}'.format(args.Abam))

            if not args.Bbam:
                self.parser.error('-Bbam is required with -w target_SNP_selection')
            elif not os.path.isfile(args.Bbam):
                self.parser.error('No such file {0}'.format(args.Bbam))

            if not args.F1bam:
                self.F1bam=""
            elif not os.path.isfile(args.F1bam):
                self.parser.error('No such file {0}'.format(args.F1bam))

            if not args.Aname:
                self.Aname="A"

            if not args.Bname:
                self.Bname="B"

            if not args.F1sim:
                self.F1sim='{0}/sim_Aa_95.txt'.format(os.path.dirname(__file__))
            elif not os.path.isfile(args.F1sim):
                self.parser.error('No such file {0}'.format(args.F1sim))

            if not args.reference:
                self.parser.error('-reference is required with -w target_SNP_selection')
            elif not os.path.isfile(args.reference):
                self.parser.error('No such file {0}'.format(args.reference))
            
            if not args.position:
                self.parser.error('-position is required with -w target_SNP_selection')
            
            if not args.output_dir:
                self.parser.error('-output_dir is required with -w target_SNP_selection')
            
            if not args.min_depth:
                self.min_depth=10
            else:
                self.min_depth=int(self.min_depth)

            if not args.max_depth:
                self.max_depth=99
            else:
                self.max_depth=int(self.max_depth)

            if not args.Bhetero:
                self.Bhetero="no"
            elif not args.Bhetero=="yes" or args.Bhetero=="no":
                self.parser.error('-Bhetero is yes or no')

            if not args.Bsim:
                self.Bsim='{0}/sim_Aa_95.txt'.format(os.path.dirname(__file__))
            elif not os.path.isfile(args.Bsim):
                self.parser.error('No such file {0}'.format(args.Bsim))

            if not args.minMQ:
                self.minMQ=0
            else:
                self.minMQ=int(self.minMQ)

            if not args.minBQ:
                self.minBQ=13
            else:
                self.minBQ=int(self.minBQ)

        elif self.work=="CAPS":

            if not args.output_dir:
                self.parser.error('-output_dir is required with -w CAPS')
            elif not os.path.isdir(args.output_dir):
                self.parser.error('No such directory {0}'.format(args.output_dir))

            if not args.restriction_enzyme:
                self.parser.error('-restriction_enzyme is required with -w CAPS')
            elif not os.path.isfile(args.restriction_enzyme):
                self.parser.error('No such file {0}'.format(args.restriction_enzyme))

            if not args.recipe:
                self.recipe='{0}/primer_recipe.txt'.format(os.path.dirname(__file__))
            elif not os.path.isfile(args.recipe):
                self.parser.error('No such file {0}'.format(args.recipe))

            if not args.PCR_max_size:
                self.PCR_max_size=1000
            else:
                self.PCR_max_size=int(self.PCR_max_size)

            if not args.PCR_min_size:
                self.PCR_min_size=500
            else:
                self.PCR_min_size=int(self.PCR_min_size)

            if not args.fragment_min_size:
                self.fragment_min_size=200
            else:
                self.fragment_min_size=int(self.fragment_min_size)

            if not args.make_html:
                self.make_html="yes"
            elif not args.make_html=="yes" or args.make_html=="no":
                self.parser.error('-make_html is yes or no')
        
        elif self.work=="ARMS_preparation":

            if not args.output_dir:
                self.parser.error('-output_dir is required with -w ARMS_preparation')
            elif not os.path.isdir(args.output_dir):
                self.parser.error('No such directory {0}'.format(args.output_dir))

            if not args.recipe:
                self.recipe='{0}/primer_recipe.txt'.format(os.path.dirname(__file__))
            elif not os.path.isfile(args.recipe):
                self.parser.error('No such file {0}'.format(args.recipe))

        elif self.work=="tetra_ARMS":

            if not args.output_dir:
                self.parser.error('-output_dir is required with -w tetra_ARMS')
            elif not os.path.isdir(args.output_dir):
                self.parser.error('No such directory {0}'.format(args.output_dir))

            if not args.recipe:
                self.recipe='{0}/primer_recipe.txt'.format(os.path.dirname(__file__))
            elif not os.path.isfile(args.recipe):
                self.parser.error('No such file {0}'.format(args.recipe))

            if not args.first_size:
                self.first_size_min=100
                self.first_size_max=500
            elif "-" in args.first_size:
                first_size_array=args.first_size.split("-")
                self.first_size_min=int(first_size_array[0])
                self.first_size_max=int(first_size_array[1])

            if not args.second_size:
                self.second_size_min=600
                self.second_size_max=1000
            elif "-" in args.second_size:
                second_size_array=args.second_size.split("-")
                self.second_size_min=int(second_size_array[0])
                self.second_size_max=int(second_size_array[1])

            if not args.make_html:
                self.make_html="yes"
            elif not args.make_html=="yes" or args.make_html=="no":
                self.parser.error('-make_html is yes or no')

        elif self.work=="tri_ARMS":

            if not args.output_dir:
                self.parser.error('-output_dir is required with -w tetra_ARMS')
            elif not os.path.isdir(args.output_dir):
                self.parser.error('No such directory {0}'.format(args.output_dir))

            if not args.recipe:
                self.recipe='{0}/primer_recipe.txt'.format(os.path.dirname(__file__))
            elif not os.path.isfile(args.recipe):
                self.parser.error('No such file {0}'.format(args.recipe))

            if not args.PCR_max_size:
                self.PCR_max_size=700
            else:
                self.PCR_max_size=int(self.PCR_max_size)

            if not args.PCR_min_size:
                self.PCR_min_size=100
            else:
                self.PCR_min_size=int(self.PCR_min_size)

            if not args.SNP_dist:
                self.SNP_dist_min=100
                self.SNP_dist_max=300
            elif "-" in args.SNP_dist:
                SNP_dist_array=args.SNP_dist.split("-")
                self.SNP_dist_min=int(SNP_dist_array[0])
                self.SNP_dist_max=int(SNP_dist_array[1])

            if not args.make_html:
                self.make_html="yes"
            elif not args.make_html=="yes" or args.make_html=="no":
                self.parser.error('-make_html is yes or no')

    def run(self):
        
        if self.work=="target_SNP_selection":
            from dnamarkmaker_script.target_SNP_selection import target_SNP_selection
            cmd=target_SNP_selection(self.Abam,self.Bbam,self.F1bam,self.Aname,self.Bname,self.F1sim,self.reference,self.position,
                                     self.output_dir,self.min_depth,self.max_depth,self.Bhetero,self.Bsim,self.minMQ,self.minBQ)
            cmd.run()
        elif self.work=="CAPS":
            from dnamarkmaker_script.CAPS import CAPS
            cmd=CAPS(self.output_dir,self.restriction_enzyme,self.recipe,self.PCR_max_size,self.PCR_min_size,self.fragment_min_size,self.thread,self.make_html)
            cmd.run()
        elif self.work=="ARMS_preparation":
            from dnamarkmaker_script.ARMS_preparation import ARMS_preparation
            cmd=ARMS_preparation(self.output_dir,self.recipe,self.thread)
            cmd.run()
        elif self.work=="tetra_ARMS":
            from dnamarkmaker_script.tetra_ARMS import tetra_ARMS
            cmd=tetra_ARMS(self.output_dir,self.recipe,self.thread,self.first_size_min,self.first_size_max,self.second_size_min,self.second_size_max,self.make_html)
            cmd.run()
        elif self.work=="tri_ARMS":
            from dnamarkmaker_script.tri_ARMS import tri_ARMS
            cmd=tri_ARMS(self.output_dir,self.recipe,self.thread,self.PCR_max_size,self.PCR_min_size,self.SNP_dist_min,self.SNP_dist_max,self.make_html)
            cmd.run()

if __name__ == '__main__':
    # tracemalloc.start()
    # print('[MarkMaker:{}] start MarkMaker'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    DNAMarkMaker().run()
    # print('[MarkMaker:{}] finish MarkMaker'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    # current, peak = tracemalloc.get_traced_memory()
    # print("Peak memory usage: {} GB".format(peak / 10**9))
    # tracemalloc.stop()

