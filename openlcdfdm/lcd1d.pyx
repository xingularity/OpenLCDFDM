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

from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "LCD_Optics2x2.hpp" namespace "LCDOptics":
    cdef enum OpticalMaterialClass:
        OPT_GLASS = 0
        OPT_ISOTROPIC = 1
        OPT_UNIAXIAL = 2
        OPT_LCMATERIAL = 3
        OPT_POLARIZER = 4
        OPT_BIAXIAL = 5

cdef extern from "LCD1D_FDM1DSolver.hpp" namespace "LCD1D":
    cdef struct LCParamters:
        double thick
        double epsr_para
        double epsr_perp
        double gamma
        double k11
        double k22
        double k33
        double q0

    cdef struct DielecParameters:
        double thick
        double epsr

    cdef struct RubbingCondition:
        double tftTheta
        double tftPhi
        double cfTheta
        double totalTwist

cdef extern from "LCD1D_LCD1DMain.hpp" namespace "LCD1D":
    cdef LCParamters createLCParameters(double _thick, double _epsr_para, double _epsr_perp, double _gamma, double _k11, double _k22, double _k33, double _q0)
    cdef DielecParameters createDielectricParameters(double _thick, double _epsr)
    cdef RubbingCondition createRubbingCondition(double _tftTheta, double _tftPhi, double _cfTheta, double _totalTwist)

    cdef cppclass LCD1DMainBase:
        void setTFTPI(DielecParameters _tftpi)
        void setCFPI(DielecParameters _cfpi)
        void setOMPThreadNum(size_t _num)
        size_t addOpticalGlassLayer(double _thick, map[double, double complex] _nkSpectrum, int pos)
        size_t addOpticalIsotropicLayer(double _thick, map[double, double complex] _nkSpectrum, OpticalMaterialClass _class, int pos)
        size_t addOpticalPolarizer(double _thick, map[double, vector[complex]] _nkSpectrum, vector[vector[double]] _axes, int pos)
        size_t addOpticalUnaixialLayer(double _thick, map[double, vector[complex]] _nkSpectrum, OpticalMaterialClass _class, vector[vector[double]] _axes , int pos)
        size_t addOpticalLC(double _thick, map[double, vector[complex]] _nkSpectrum, int pos)
        void removeOpticalLayer(size_t _index)
        void setOpticalIncidentAngles()
        void setOpticalIncidentAngles(double _thetaInterval, double _phiInterval)
        void setOpticalIncidentAngles(vector[pair[double, double]] _angles)
        void setOpticalWavelength(double _lambda_start, double _lambda_end, double _lambda_step)
        void setOpticalWavelength(double _lambda)
        void setOpticalSourceSpectrum(map[double, double] _input)
        void resetLCParam(const LCParamters _param, const size_t _lcLayerNum, double _dt)
        void resetLCRubbing(const RubbingCondition _rubbing);
        void useOptical2X2Lambertian(bool _if);
        void createExtendedJones();
        vector[vector[pair[double, double]]] getIncidentAngles()const;

    cdef cppclass LCD1DStaticMain(LCD1DMainBase):
        LCD1DStaticMain(double _lcLayerNum, double _dt, LCParamters _lcParam, RubbingCondition _rubbing, \
            double _voltStart, double _voltEnd, double _voltStep, double _maxIter, double _maxError);
        LCD1DStaticMain()
        vector[double] getCalcVolts() const
        vector[vector[vector[double]]] getTransmissions()const
        vector[vector[vector[double]]] getLCDirResults()const
        vector[double] getNormalTransmissions()const
        void resetCalcVolts(double _voltStart, double _voltEnd, double _voltStep)
        void calculate()

cdef class pyLCD1DStaticMain:
    cdef LCD1DStaticMain *thisptr
    cdef LCParamters lcParamStruct
    cdef RubbingCondition rubbingCondSruct
    cdef DielecParameters dielecParamStruct

    def __cinit__(self, **kwargs):
        if (len(kwargs) == 0):
            self.thisptr = new LCD1DStaticMain()
        else:
            lcLayerNum = kwargs['lcLayerNum']
            dt = kwargs['dt']
            voltStart = kwargs['voltStart']
            voltEnd = kwargs['voltEnd']
            voltStep = kwargs['voltStep']
            maxIter = kwargs['maxIter']
            maxError = kwargs['convergeError']
            lcparam = kwargs['lcparam']
            rubbing = kwargs['rubbing']
            lcParamStruct = createLCParameters(lcparam['thick'], lcparam['epsr_para'], lcparam['epsr_perp'], lcparam['gamma'], lcparam['k11'],lcparam['k22'],lcparam['k33'],lcparam['q0'])
            rubbingCondSruct = createRubbingCondition(rubbing['tftTheta'], rubbing['tftPhi'], rubbing['cfTheta'], rubbing['totalTwist'])
            self.thisptr = new LCD1DStaticMain(lcLayerNum, dt, lcParamStruct, rubbingCondSruct, voltStart, voltEnd, voltStep, maxIter, maxError)

    def __dealloc__(self):
        del self.thisptr

    def setTFTPI(self, **kwargs):
        if (len(kwargs) == 0):
            dielecParamStruct = createDielectricParameters(0.0, 0.0)
        else:
            dielecParamStruct = createDielectricParameters(kwargs['thick'], kwargs['epsr'])
        self.thisptr.setTFTPI(dielecParamStruct)

    def setOMPThreadNum(self, size_t _num):
        self.thisptr.setOMPThreadNum(_num)

    def addOpticalGlassLayer(self, _thick not None, dict _nkSpectrum not None, int pos=None):
        if (pos == None):
            pos = -1
        return self.thisptr.addOpticalGlassLayer(_thick, _nkSpectrum, pos)

    def addOpticalIsotropicLayer(self, _thick, dict _nkSpectrum not None, _class, int pos):
        if (pos == None):
            pos = -1
        return self.thisptr.addOpticalIsotropicLayer(_thick, _nkSpectrum, _class, pos)

    def addOpticalPolarizer(self, _thick, dict _nkSpectrum not None, _axes, int pos):
        if (pos == None):
            pos = -1
        return self.thisptr.addOpticalPolarizer(_thick, _nkSpectrum, _axes, pos)

    def addOpticalUnaixialLayer(self, _thick, dict _nkSpectrum not None, _class, _axes , int pos):
        if (pos == None):
            pos = -1
        return self.thisptr.addOpticalUnaixialLayer(_thick, _nkSpectrum, _class, _axes, pos)

    def addOpticalLC(self, _thick, dict _nkSpectrum not None, int pos):
        if (pos == None):
            pos = -1
        return self.thisptr.addOpticalLC(_thick, _nkSpectrum, pos)

    def removeOpticalLayer(self, _index):
        self.thisptr.removeOpticalLayer(_index)

    def setOpticalIncidentAnglesWithInterval(self, _thetaInterval, _phiInterval):
        self.thisptr.setOpticalIncidentAngles(_thetaInterval, _phiInterval)
