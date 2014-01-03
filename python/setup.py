# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 13:50:22 2013

@author: jbalzer
"""

from distutils.core import setup, Extension
import numpy as np

# define c-extensions
r4rio = Extension('r4r.core._io',
                         include_dirs = ['/usr/include/',np.get_include(),'../r4r_core'],
                         sources=['r4r/core/io.cpp'],
                         library_dirs=['../../build/r4r_core'],
                         libraries=['r4r_core'],
                         extra_compile_args=['-std=c++0x'])

r4rsplines = Extension('r4r.core._bsplines',
                         include_dirs = ['/usr/include/',np.get_include(),'../r4r_core'],
                         sources=['r4r/core/bsplines.cpp'],
                         library_dirs=['../../build/r4r_core'],
                         libraries=['r4r_core'],
                         extra_compile_args=['-std=c++0x'])


# the core setup routine
setup(name='r4r',
      version='1.0',
      description='R4R Python interface.',
      author='Jonathan Balzer',
      author_email='jonabalzer@gmail.com',
      url='https://sites.google.com/site/jonabalzer/',
      packages=['r4r','r4r.core'],
      ext_modules=[r4rio,r4rsplines])

# MAKE SURE NOT TO RUN PYTHON INTERPRETER FROM THE DIRECTORY THAT HOLDS SETUP.
# THIS WILL MAKE THE C-EXTENSIONS INVISIBLE!!!