#########################################################################
## Copyright (C) 2015 Zong-han, Xie <icbm0926@gmail.com>.
## All rights reserved.
##
## You may use this file under the terms of the BSD license as follows:
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##   * Neither the name of OpenLCDFDM nor the names of its contributors
##     may be used to endorse or promote products derived from this
##     software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
##
#########################################################################
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

sourcefiles=['./openlcdfdm/lcd1d.pyx', './src/LCD1D_ExtendedJones.cpp', \
    './src/LCD1D_FDM1DSolver.cpp', './src/LCD1D_LCD1DMain.cpp', './src/LCD_Optics2x2.cpp', \
    './src/LCD_SpectrumInterpolator.cpp', './src/LCD_TimeWaveform.cpp']

import platform as pm
import sys
if pm.system() == 'Darwin':
    sitepkg_path=""
    for i in sys.path:
        if i.count("site-packages")>0:
            sitepkg_path=i
    extensions = [Extension('lcd', sourcefiles, language="c++", extra_compile_args=['-std=c++11', '-O3', '-fopenmp'], include_dirs = [sitepkg_path+"/numpy/core/include/", "./include/"])]
else:
	extensions = [Extension('lcd', sourcefiles, language="c++", extra_compile_args=['-std=c++11', '-O3', '-fopenmp'], include_dirs = ["./include/"])]

setup(ext_modules = cythonize(extensions)
      )
