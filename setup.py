#!/usr/bin/env python3

from distutils.core import setup
from dnamarkmaker_script.__init__ import __version__

setup(name='DNAMarkMaker',
      version='{}'.format(__version__),
      description='DNAMarkMaker: pipline to develop ARMS and CAPS marker',
      author_name='Tenta Segawa',
      url='https://github.com/SegawaTenta/DNAMarkMaker-CUI',
      packages=['dnamarkmaker_script'],
      install_requires=["primer3-py"],
      entry_points={'console_scripts':['DNAMarkMaker=dnamarkmaker_script.DNAMarkMaker:DNAMarkMaker']}
      )