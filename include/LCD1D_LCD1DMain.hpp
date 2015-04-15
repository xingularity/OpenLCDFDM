/*
 * Copyright (C) 2015 Zong-han, Xie <icbm0926@gmail.com>.
 *
 * You may use this file under the terms of the BSD license as follows:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of OpenLCDFDM nor the names of its contributors
 *     may be used to endorse or promote products derived from this
 *     software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
 *
 */

#ifndef LCD1D_LCD1DMAIN_HPP
#define LCD1D_LCD1DMAIN_HPP

#include <omp.h>
#include <ctime>
#include <stdexcept>
#include "LCD1D_ExtendedJones.hpp"
#include "LCD_UsefulFuncs.hpp"
#include "LCD1D_FDM1DSolver.hpp"

class LCD1DMainBase{
public:
    void setTFTPI(LCD1D::DielecParameters _tftpi);
    void setCFPI(LCD1D::DielecParameters _cfpi);
    void setOMPThreadNum(size_t _num);
    void addOpticalGlassLayer(double _thick, std::vector<double> _spectrumLambdas, std::vector<std::complex<double> > _nkSpectrum);
    void addOpticalIsotropicLayer(double _thick, std::vector<double> _spectrumLambdas, 
        std::vector<std::complex<double> > _nkSpectrum, OpticalMaterialClass _class);
    void addOpticalPolarizer(double _thick, std::vector<double> _spectrumLambdas, std::vector<std::complex<double> > _nokoSpectrum, 
        std::vector<std::complex<double> > _nekeSpectrum);
    void addOpticalUnaixialLayer(double _thick, std::vector<double> _spectrumLambdas, 
        std::vector<std::complex<double> > _nokoSpectrum, std::vector<std::complex<double> > _nekeSpectrum, 
        OpticalMaterialClass _class);
    void addOpticalLC(double _thick, std::vector<double> _spectrumLambdas, std::vector<std::complex<double> > _nokoSpectrum, 
        std::vector<std::complex<double> > _nekeSpectrum);
    ///if input 0, then only calculate normal incident.
    void setOpticalIncidentAngles(size_t _intervalDegree);
    void setOpticalIncidentAngles(std::vector<std::pair<double, double> > _angles);
    void setOpticalLambda(std::vector<double> _lambdas);
    void setOpticalLambda(double _lambda);
    void setOpticalSourceSpectrum(std::vector<double> _lambdas, std::vector<double> _powers);
    void calculate();
protected:
    LCD1DMainBase(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing);
    std::shared_ptr<LCD1D::Epsilon> epsilonr;
    std::shared_ptr<LCD1D::Potential> potentials;
    std::shared_ptr<LCD1D::LCDirector> lcDir;
    LCDOptics::MATERIALLAYERS2X2CONT materials;
    std::shared_ptr<LCDOptics::ExtendedJones> extj;
    LCDOptics::IAngles inAngles;
    LCD::DOUBLEARRAY1D lambdas;
};

/**
Doing static LCD simulation
*/
class LCD1DStaticMain: public LCD1DMainBase{
public:
    LCD1DStaticMain(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing, 
        double _voltStart, double _voltEnd, double _voltStep, double _maxIter, double _error);
    ///vector records which voltages are calculated.
    LCD::DOUBLEARRAY1D getCalcVolts() const;
    ///std::vector 3D container. [voltage index][iAngle index][iAngle index]
    std::vector<LCDOptics::TRANSRESULT> getTransmissions()const;
    ///incident angles. originally using tuples, return std::pairs for Cython.
    std::vector<std::vector<std::pair<double, double> > > getIncidentAngles()const;
    ///[volts index][z-grid index][component index]
    std::vector<DOUBLEARRAY2D> getLCDirResults()const;

private:
    double maxIter;
    double error;
    LCD::DOUBLEARRAY1D  calcVolts;
    std::vector<DOUBLEARRAY2D> lcDirResult;
    std::vector<LCDOptics::TRANSRESULT> transResults;
};

/**
Doing dynamic LCD simulation
*/
class LCD1DDynamicMain: public LCD1DMainBase{
public:
    LCD1DStaticMain(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing, double _maxCalcTime);
    ///set DC wavreform if one want to use DC voltage B.C. to simulate.
    void setDCWaveform(double _volt);
    ///set step waveform if one want to use step waveform profile to simulate.
    void setStepWaveform(std::map<double, double> _profile);
    ///vector records which step data are recorded at.
    LCD::DOUBLEARRAY1D getRecordStep() const;
    ///vector records which time (in ms) data are recorded at.
    LCD::DOUBLEARRAY1D getRecordTime() const;
    ///std::vector 3D container. [record step index][iAngle index][iAngle index]
    std::vector<LCDOptics::TRANSRESULT> getTransmissions()const;
    ///incident angles. originally using tuples, return std::pairs for Cython.
    std::vector<std::vector<std::pair<double, double> > > getIncidentAngles()const;
    ///[record time index][z-grid index][component index]
    std::vector<DOUBLEARRAY2D> getLCDirResults()const;
private:
    std::shared_ptr<WaveformBase> waveform;
    double maxCalcTime;
    LCD::DOUBLEARRAY1D recordSteps;
    LCD::DOUBLEARRAY1D recordTimes;
    std::vector<DOUBLEARRAY2D> lcDirResult;
    std::vector<LCDOptics::TRANSRESULT> transResults;
};

#endif
