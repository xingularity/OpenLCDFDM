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

    class LCDirector{
        friend class LCSolver;
    public:
        ///Constructor
        LCDirector(const LCParamters lcParam_, const RubbingCondition rubbing_, size_t layerNum);
        size_t getSize()const{return dir.size();}
        const DIRVEC getDirectors()const {return lcDir;}
        LCParamters getLCParams()const {return lcParam;}
        RubbingCondition getLCRubbing()const {return rubbing;}
        double getDirectors(int iz, int i)const {return lcDir(iz)(i);}
        double& getDirectors(int iz, int i){return lcDir(iz)(i);}
        BzVECD3D& getDirectors(int iz){return lcDir(iz);}
        const BzVECD3D getDirectors(int iz)const {return lcDir(iz);}
        const DOUBLEARR1D& getEps33()const {return eps33;}
        ///reset LC directors according to the rubbing BC.
        void resetLCDirectors();
        ///reset LC parameters
        void resetConditions(const LCParamters lcParam_);
        ///reset rubbing conditions
        void resetConditions(const RubbingCondition rubbing_);
    private:
        ///LC parameters
        LCParamters lcParam;
        ///rubbing conditions
        RubbingCondition rubbing;
        ///LC directors
        LCD::DIRVEC lcDir;
        ///
        DOUBLEARR1D eps33;
        std::shared_ptr<LCSolver> solver;
    };

    ///only for potentials in LC, potentials between pi and LC will be calculated with capacitance.
    class Potential{
        friend class PotentialSolver1D;
    public:

        Potential(size_t size_);
        const getSize();
        const BZARRAY<double, 1> getPotentials()const ;
        BZARRAY<double, 1> getPotentials();
        double& getPotentials(size_t i);
        const double& getPotentials(size_t i)const ;
    private:
        BZARRAY<double, 1> potential;
    };

    class SolverBase{
    public:
        SolverBase(){};
        SolverBase(double _dz):dz(_dz){};
        ///move one step forward.
        virtual void update(double t) = 0;
    protected:
        double dz{0.0};
    };

    class PotentialSolver1D:public SolverBase{
    public:
        PotentialSolver1D(Potnetial& _pot, LCDirector& _lcDir, double _dz, DielecParameters _tftpi, DielecParameters _cfpi,
            std::shared_ptr<WaveformBase> voltWavePtr_ = std::shared_ptr<LCD::WaveformBase>());
        virtual void update(double t);
    protected:
        DielecParameters tftpi;
        DielecParameters cfpi;
        Potential& potentials;
        LCDirector& lcDir;
        Epsilon& epsilon;
        std::shared_ptr<LCD::WaveformBase> voltWavePtr;
    };

    class LCSovler:public SolverBase{
    public:
        LCSolver(Potnetial& pot_, LCDirector& lcDir_, LCParamters _lcParam, double _dz, double _dt);
        virtual void update(double t);
    protected:
        ///LC parameters
        LCParamters lcParam;
        Potential& potentials;
        LCDirector& lcDir;
    };
};

#endif
