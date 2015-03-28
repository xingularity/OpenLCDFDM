/*
 * Copyright (C) 2015 Zong-han, Xie <icbm0926@gmail.com>.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef LCD1D_FDM1DSOLVER_HPP
#define LCD1D_FDM1DSOLVER_HPP

#include "LCD_ContainerDefine.hpp"
#include "LCD_TimeWaveform.hpp"
#include <memory>

namespace LCD1D{
    class SolverBase;
    class PotentialSolver1D;
    class LCSolver;
    class TimeWaveform;

    ///Rubbing conditions are tft-side theta, tft-side phi, cf-side theta, total twist.
    struct LCParamters{
        double thick;
        double epsPara;
        double epsPerp;
        double gamma;
        double k11;
        double k22;
        double k33;
    };

    struct DielecParameters{
        double thick;
        double eps;
    };

    ///Rubbing condition
    struct RubbingCondition{
        double tftTheta;
        double tftPhi;
        double cfTheta;
        double totalTwist;
    };

    /**
    This class only records the epsilon distribution in LC.
    Epsilon values of tftpi & cfpi are recorded with the values only.
    */
    class Epsilon{
    public:
        Epsilon(size_t _lc_layernum){lcLayerNum=_lc_layernum;}
        ///set TFTPI parameters
        void setTFTPI(DielecParameters _param){tftpi_eps = _param.eps; tftpiThick = _param.thick;}
        ///set CFPI parameters
        void setCFPI(DielecParameters _param){cfpi_eps = _param.eps; cfpiThick = _param.thick;}
        ///get epsilon distribution in LC
        DOUBLEARRAY1D getLCEpsilon()const{return eps33;}
        ///get epsilon distribution in LC
        DOUBLEARRAY1D& getLCEpsilon(){return eps33;}
        ///get epsilon with an index in LC
        double& getLCEpsilon(size_t i){return eps33[i];}
        ///get epsilon with an index in LC
        double getLCEpsilon(size_t i)const{return eps33[i];}
        /// calculate cell capacitance, units in pF/cm^2
        double calculateCapacitance(double dz) const ;
        ///get TFTPI capacitance, units in pF/cm^2
        double getTFTCapcitance()const {return tftCap;}
        ///get CFPI capacitance, units in pF/cm^2
        double getCFCapcitance()const {return cfCap;}

    private:
        size_t lcLayerNum{0};
        DOUBLEARRAY1D eps33;
        double tftpi_eps{0.0};
        double cfpi_eps{0.0};
        double tftpiThick{0.0};
        double cfpiThick{0.0};
        ///capacitance of tftpi
        double tftCap{0.0};
        ///capacitance of cfpi
        double cfCap{0.0};
    };

    /**
    record LC directors
    */
    template <class UpdatePolicy>
    class LCDirector{
        ///C++11 syntax allowed
        friend UpdatePolicy;
    public:
        ///Constructor
        LCDirector(const LCParamters _lcParam, const RubbingCondition _rubbing, Epsilon& _epsilion);
        size_t getSize()const{return dir.size();}
        const DIRVEC getDirectors()const {return lcDir;}
        LCParamters getLCParams()const {return lcParam;}
        RubbingCondition getLCRubbing()const {return rubbing;}
        double getDirectors(int iz, int i)const {return lcDir(iz)(i);}
        const BzVECD3D getDirectors(int iz)const {return lcDir(iz);}
        const DOUBLEARR1D& getEps33()const {return eps33;}
        ///reset LC directors according to the rubbing BC.
        void resetLCDirectors();
        ///reset LC parameters
        void resetConditions(const LCParamters _lcParam);
        ///reset rubbing conditions
        void resetConditions(const RubbingCondition _rubbing);
        ///
    private:
        ///LC parameters
        LCParamters lcParam;
        ///rubbing conditions
        RubbingCondition rubbing;
        ///LC directors
        LCD::DIRVEC lcDir;
        ///epsilon values
        Epsilon& epsilon;
        ///The functor used to update LC directors
        std::shared_ptr<UpdatePolicy> lcUpdater;
    };

    ///only for potentials in LC, potentials between pi and LC will be calculated with capacitance.
    template <class UpdatePolicy>
    class Potential{
        friend class UpdatePolicy;
    public:
        Potential(size_t size_, const Epsilon& epsilon);
        const getSize();
        const DOUBLEARRAY1D getPotentials()const ;
        DOUBLEARRAY1D getPotentials();
        const double& getPotentials(size_t i)const;
        void update(double t);
    private:
        ///epsilons
        const Epsilon& epsilons;
        /// potnetials inside LC
        DOUBLEARRAY1D potential;
        ///The functor used to update potentials
        std::shared_ptr<UpdatePolicy> potentialUpdater;
    };

    class PotentialSolversForStatic{
    public:
        PotentialSolver1D(Potnetial& _pot, LCDirector& _lcDir, double _dz, DielecParameters _tftpi, DielecParameters _cfpi);
        void changeBCVolt(double _dc);
        void update(double t);
    protected:
        double dz{0.0};
        Potential& potentials;
        Epsilon& epsilon;
    };

    class PotentialSolversForDynamic{
    public:
        PotentialSolver1D(Potnetial& _pot, LCDirector& _lcDir, double _dz, DielecParameters _tftpi,
            std::shared_ptr<LCD::WaveformBase> _voltWavePtr = std::shared_ptr<LCD::WaveformBase>());
        void update(double t);
    protected:
        double dz{0.0};
        Potential& potentials;
        Epsilon& epsilon;
        std::shared_ptr<LCD::WaveformBase> voltWavePtr;
    };

    class LCSolverBase{
    public:
        void changeLCParams(LCParamters _param){lcParam = _param;}
    protected:
        LCSolverBase(Potnetial& _pot, LCDirector& _lcDir, LCParamters _lcParam, double _dz, double _dt);
        double dz;
        double dt;
        DOUBLEARR1D temp;
        ///LC parameters
        LCParamters lcParam;
        Potential& potentials;
        LCDirector& lcDir;
    }

    class LCVecUpdate:public LCSolverBase{
    public:
        LCSolver(Potnetial& _pot, LCDirector& _lcDir, LCParamters _lcParam, double _dz, double _dt):
            LCSolverBase(_pot, _lcDir, _lcParam, _dz, _dt){}
        void update(double t);
        double residuals();
    };
};

#endif
