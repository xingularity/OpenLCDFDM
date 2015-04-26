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
    ///set tftpi parameters
    void setTFTPI(LCD1D::DielecParameters _tftpi={0.0, 0.0});
    ///set cfpi parameters
    void setCFPI(LCD1D::DielecParameters _cfpi={0.0,0.0});
    ///set openMP thread numbers
    void setOMPThreadNum(size_t _num);
    void addOpticalGlassLayer(double _thick, std::map<double, std::complex<double> > _nkSpectrum);
    void addOpticalIsotropicLayer(double _thick, std::map<double, std::complex<double> > _nkSpectrum, LCDOptics::OpticalMaterialClass _class);
    void addOpticalPolarizer(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum);
    void addOpticalUnaixialLayer(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum,
        LCDOptics::OpticalMaterialClass _class);
    void addOpticalLC(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum);
    ///set incident angle in the optical calculation to normal incident only
    void setOpticalIncidentAngles();
    ///Using scan interval degree to decide which angles to do.
    void setOpticalIncidentAngles(unsigned double _thetaInterval, unsigned double _phiInterval);
    ///input multiple incident angles without fix intervlas
    void setOpticalIncidentAngles(std::vector<std::pair<double, double> > _angles);
    ///for multiwavelength calculation
    void setOpticalWavelength(double _lambda_start, double _lambda_end, double _lambda_step);
    ///for single wabelength calculation
    void setOpticalWavelength(double _lambda);
    void setOpticalSourceSpectrum(LCDOptics::LIGHTSPECTRUMDATA _input);
    void enableOptical2X2Calculation(bool _ifDo=true);
    virtual void calculate() = 0;
protected:
    ///This constructor will calculate LC
    LCD1DMainBase(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing);
    ///This constructor will not calculate LC (No FDM calculation).
    LCD1DMainBase();
    void createExtendedJones();
    std::shared_ptr<LCD1D::Epsilon> epsilonr;
    std::shared_ptr<LCD1D::Potential> potentials;
    std::shared_ptr<LCD1D::LCDirector> lcDir;
    std::shared_ptr<LCDOptics::ExtendedJones> extJonesMain;
    LCDOptics::MATERIALLAYERS2X2CONT materials;
    LCDOptics::IAngles inAngles;
    LCDOptics::LIGHTSPECTRUMDATA lightSrcSpectrum;
    ///They are lambda_start, lambda_end, lambda_step, if lambda_step == 0, it's single wavelength calculation.
    ///If all three are zero, no optical calculation should be performed.
    std::tuple<double, double, double> multiWavelengthLambdas = {0.0,0.0,0.0};
    bool ifCalculate2X2Optics = false;
};

/**
Doing static LCD simulation
*/
class LCD1DStaticMain: public LCD1DMainBase{
public:
    ///with LC calculation
    LCD1DStaticMain(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing,
        double _voltStart, double _voltEnd, double _voltStep, double _maxIter, double _error);
    ///No LC calculation, optical calculation only
    LCD1DStaticMain();
    ///vector records which voltages are calculated.
    LCD::DOUBLEARRAY1D getCalcVolts() const;
    ///std::vector 3D container. [voltage index][iAngle index][iAngle index]
    std::vector<LCDOptics::TRANSRESULT> getTransmissions()const;
    ///incident angles. originally using tuples, return std::pairs for Cython.
    std::vector<std::vector<std::pair<double, double> > > getIncidentAngles()const;
    ///[volts index][z-grid index][component index]
    std::vector<DOUBLEARRAY2D> getLCDirResults()const;
    virtual void calculate();
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
    virtual void calculate();
private:
    std::shared_ptr<WaveformBase> waveform;
    double maxCalcTime;
    LCD::DOUBLEARRAY1D recordSteps;
    LCD::DOUBLEARRAY1D recordTimes;
    std::vector<DOUBLEARRAY2D> lcDirResult;
    std::vector<LCDOptics::TRANSRESULT> transResults;
};

#endif
