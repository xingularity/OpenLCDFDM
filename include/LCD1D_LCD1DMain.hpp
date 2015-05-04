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
#include <set>
#include "LCD1D_ExtendedJones.hpp"
#include "LCD_UsefulFuncs.hpp"
#include "LCD1D_FDM1DSolver.hpp"

namespace LCD1D{

    class LCD1DMainBase{
    public:
        ///set tftpi parameters
        void setTFTPI(LCD1D::DielecParameters _tftpi={0.0, 0.0});
        ///set cfpi parameters
        void setCFPI(LCD1D::DielecParameters _cfpi={0.0,0.0});
        ///set openMP thread numbers
        void setOMPThreadNum(size_t _num);
        size_t addOpticalGlassLayer(double _thick, std::map<double, std::complex<double> > _nkSpectrum, int pos = -1);
        size_t addOpticalIsotropicLayer(double _thick, std::map<double, std::complex<double> > _nkSpectrum, LCDOptics::OpticalMaterialClass _class, int pos = -1);
        size_t addOpticalPolarizer(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum, LCD::DOUBLEARRAY2D _axes, int pos = -1);
        size_t addOpticalUnaixialLayer(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum,
            LCDOptics::OpticalMaterialClass _class, LCD::DOUBLEARRAY2D _axes = LCD::DOUBLEARRAY2D(), int pos = -1);
        size_t addOpticalLC(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum, int pos = -1);
        void removeOpticalLayer(size_t _index);
        ///set incident angle in the optical calculation to normal incident only
        void setOpticalIncidentAngles();
        ///Using scan interval degree to decide which angles to do.
        void setOpticalIncidentAngles(double _thetaInterval, double _phiInterval);
        ///input multiple incident angles without fix intervlas
        void setOpticalIncidentAngles(std::vector<std::pair<double, double> > _angles);
        ///for multiwavelength calculation
        void setOpticalWavelength(double _lambda_start, double _lambda_end, double _lambda_step);
        ///for single wabelength calculation
        void setOpticalWavelength(double _lambda);
        void setOpticalSourceSpectrum(LCDOptics::LIGHTSPECTRUMDATA _input);
        void resetLCParam(const LCD1D::LCParamters _param, const size_t _lcLayerNum, double _dt=0.0);
        void resetLCRubbing(const LCD1D::RubbingCondition _rubbing);
        void useOptical2X2Lambertian(bool _if=true);
        void createExtendedJones();
        ///incident angles. originally using tuples, return std::pairs for Cython.
        std::vector<std::vector<std::pair<double, double> > > getIncidentAngles()const;
        virtual void calculate() = 0;
    protected:
        ///This constructor will calculate LC
        LCD1DMainBase(double _lcLayerNum, double _dt, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing);
        ///This constructor will not calculate LC (No FDM calculation).
        LCD1DMainBase();
        std::shared_ptr<LCD1D::Epsilon> epsilonr;
        std::shared_ptr<LCD1D::Potential> potentials;
        std::shared_ptr<LCD1D::LCDirector> lcDir;
        std::shared_ptr<LCD1D::ExtendedJones> extJonesMain;
        LCDOptics::MATERIALLAYERS2X2CONT materials;
        LCDOptics::IAngles inAngles;
        LCDOptics::LIGHTSPECTRUMDATA lightSrcSpectrum;
        ///They are lambda_start, lambda_end, lambda_step, if lambda_start == lambda_end, it's single wavelength calculation.
        ///default calculate 550nm single wavelength.
        std::tuple<double, double, double> multiWavelengthLambdas = std::make_tuple(0.55,0.55,0.0);
        bool ifCalculate2X2Optics = false;
        bool ifUseLambertian = false;
        double dt;
    };

    /**
    Doing static LCD simulation
    */
    class LCD1DStaticMain: public LCD1DMainBase{
    public:
        ///with LC calculation
        LCD1DStaticMain(double _lcLayerNum, double _dt, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing,
            double _voltStart, double _voltEnd, double _voltStep, double _maxIter, double _maxError);
        ///No LC calculation, optical calculation only
        LCD1DStaticMain();
        ///vector records which voltages are calculated.
        LCD::DOUBLEARRAY1D getCalcVolts() const;
        ///std::vector 3D container. [voltage index][iAngle index][iAngle index]
        std::vector<LCD1D::TRANSRESULT> getTransmissions()const;
        ///[volts index][z-grid index][component index]
        std::vector<LCD::DOUBLEARRAY2D> getLCDirResults()const;
        //return the records of the transmissions of normal incident lights.
        LCD::DOUBLEARRAY1D getNormalTransmissions()const;
        ///reset scanning voltages
        void resetCalcVolts(double _voltStart, double _voltEnd, double _voltStep);
        ///main function
        virtual void calculate();

    private:
        double calculateOneVolt(double _volt);
        void calc2X2OpticsOneSetLCDir(LCD::DIRVEC lcDir);
        double maxIter;
        double maxError;
        LCD::DOUBLEARRAY1D calcVolts;
        std::vector<LCD::DOUBLEARRAY2D> lcDirResults;
        std::vector<LCD1D::TRANSRESULT> transResults;
        LCD::DOUBLEARRAY1D normalTransmissions;
    };

    /**
    Doing dynamic LCD simulation
    */
    class LCD1DDynamicMain: public LCD1DMainBase{
    public:
        LCD1DDynamicMain(double _lcLayerNum, double _dt, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing, double _maxCalcTime);
        ///Directors and transmissions will be calculated and recorded at these time steps
        void setRecordTime(LCD::DOUBLEARRAY1D steps);
        ///Directors and transmissions will be calculated and recorded at these time steps
        void setRecordInterval(double _interval);
        ///set DC wavreform if one want to use DC voltage B.C. to simulate.
        void setDCWaveform(double _volt);
        ///set step waveform if one want to use step waveform profile to simulate.
        void setStepWaveform(std::map<double, double> _profile, double period);
        ///vector records which step data are recorded at.
        std::vector<size_t> getRecordStep() const;
        ///vector records which time (in ms) data are recorded at.
        LCD::DOUBLEARRAY1D getRecordTime() const;
        ///std::vector 3D container. [records index][iAngle index][iAngle index]
        std::vector<LCD1D::TRANSRESULT> getTransmissions()const;
        ///[record time index][z-grid index][component index]
        std::vector<LCD::DOUBLEARRAY2D> getLCDirResults()const;
        //return the records of the transmissions of normal incident lights.
        LCD::DOUBLEARRAY1D getNormalTransmissions()const;
        virtual void calculate();
    private:
        ///record steps
        std::set<unsigned long> rsteps;
        void checkStepsToDumpAndCalcOptics(size_t iternum);
        void calc2X2OpticsOneSetLCDir(LCD::DIRVEC directors);
        std::shared_ptr<LCD::WaveformBase> waveform;
        double maxCalcTime;
        std::vector<size_t> recordSteps;
        LCD::DOUBLEARRAY1D recordTimes;
        std::vector<LCD::DOUBLEARRAY2D> lcDirResults;
        std::vector<LCD1D::TRANSRESULT> transResults;
        LCD::DOUBLEARRAY1D normalTransmissions;
    };
};
#endif
