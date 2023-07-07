#!/usr/bin/env python3

from distutils.core import setup
from markmaker_script.__init__ import __version__

setup(name='MarkMaker',
      version='{}'.format(__version__),
      description='MarkMaker: pipline to develop ARMS and CAPS marker',
      author_name='Tenta Segawa',
      url='https://github.com/SegawaTenta/MarkMaker-CUI',
      package=['MarkMaker'],
      entry=points={'console_scripts':['MarkMaker=markmaker_script.MarkMaker:main']}
      )