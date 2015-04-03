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

#include "LCD_ConstantDefine.hpp"
#include "LCD_ContainerDefine.hpp"
#include "LCD_TimeWaveform.hpp"
#include <memory>
#include <stdexcept>

namespace LCD1D{
    using LCD::DIRVEC;
    using LCD::DOUBLEARRAY1D;
    using LCD::BzVECD3D;

    class LCSolverBase;
    class Epsilon;
    class Potential;
    class LCDirector;
    class PotentialCalculate;
    class PotentialSolversForStatic;
    class PotentialSolversForDynamic;
    class LCVecUpdater;

    ///Rubbing conditions are tft-side theta, tft-side phi, cf-side theta, total twist.
    struct LCParamters{
        double thick;
        double epsr_para;
        double epsr_perp;
        double gamma;
        double k11;
        double k22;
        double k33;
        double q0;
    };

    struct DielecParameters{
        double thick;
        double epsr;
    };

    ///Rubbing condition
    struct RubbingCondition{
        double tftTheta;
        double tftPhi;
        double cfTheta;
        double totalTwist;
    };

    /**
    This class only records the epsilon distribution in LC. <br/>
    Epsilon values of tftpi & cfpi are recorded with the values only. <br/>
    Note: <br/>
    EPS0 is in the unit of pF/m^2. Capacitance per unit area is calculated to be in the unit of
    1.0e-10 F/cm^2 with dz given in the unit of um.
    All the record capacitance values in this class are in the unit of 1.0e-10F/cm^2. <br/>
    */
    class Epsilon{
    public:
        Epsilon(size_t _lc_layernum, double _cellgap){
            lcLayerNum=_lc_layernum; epsr33.resize(lcLayerNum); dz=_cellgap/_lc_layernum;
        }
        ///set TFTPI parameters
        void setTFTPI(DielecParameters _param){tftpi_epsr = _param.epsr; tftpiThick = _param.thick;}
        ///set CFPI parameters
        void setCFPI(DielecParameters _param){cfpi_epsr = _param.epsr; cfpiThick = _param.thick;}
        ///get epsilon distribution in LC
        const DOUBLEARRAY1D& getLCEpsr33()const{return epsr33;}
        ///update epsr_33 in LC
        void updateEpsilonr(const double& epsr_para, const double& epsr_perp, const DIRVEC dirs);
        ///get epsilon distribution in LC
        DOUBLEARRAY1D& getLCEpsr33(){return epsr33;}
        /// calculate cell capacitance, in the unit of 1.0e-10F/cm^2
        double getLCLayerCapacitance()const{return lcCap;}
        ///get TFTPI capacitance, the unit is 1.0e-10F/cm^2
        double getTFTCapcitance()const {return tftpiThick/tftpi_epsr;}
        ///get CFPI capacitance, the unit is 1.0e-10F/cm^2
        double getCFCapcitance()const {return cfpiThick/cfpi_epsr;}
        ///get the voltage B.C. for the LC layer.
        double getLCVoltRatio()const;
    private:
        /// calculate LC cell capacitance, in the unit of 1.0e-10F/cm^2
        double calculateCapacitance() ;
        size_t lcLayerNum{0};
        double dz{0.0};
        ///eps33 in the LC layer.
        DOUBLEARRAY1D epsr33;
        double tftpi_epsr{0.0};
        double cfpi_epsr{0.0};
        double tftpiThick{0.0};
        double cfpiThick{0.0};
        ///capacitance for LC layer
        double lcCap{0.0};
        ///get total capacitance
        double totalCap{0.0};
    };

    /**
    record LC directors
    */
    class LCDirector{
        friend LCVecUpdater;
    public:
        ///Constructor
        LCDirector(size_t _lcLayerNum, double _cellgap, const LCParamters _lcParam, const RubbingCondition _rubbing, Epsilon& _epsilonr)
            :layerNum(_lcLayerNum), cellgap(_cellgap), lcParam(_lcParam), rubbing(_rubbing), epsilonr(_epsilonr){
            
            lcDir.resizeAndPreserve(_lcLayerNum+1);
            this->resetLCDirectors();
        }
        size_t getSize()const{return lcDir.extent(0);}
        size_t getLayerNum()const{return lcDir.extent(0) - 1;}
        const DIRVEC getDirectors()const {return lcDir;}
        LCParamters getLCParams()const {return lcParam;}
        RubbingCondition getLCRubbing()const {return rubbing;}
        double getDirectors(int iz, int i)const {return lcDir(iz)(i);}
        const BzVECD3D getDirectors(int iz)const {return lcDir(iz);}
        ///reset LC directors according to the rubbing BC.
        void resetLCDirectors(){
            size_t layerNum = lcDir.extent(0)-1;
            double dtheta = (rubbing.cfTheta - rubbing.tftTheta)/layerNum;
            double dphi = rubbing.totalTwist/layerNum;
            for (int i = 0; i < layerNum + 1; ++i){
                lcDir(i)(0) = std::sin(rubbing.tftTheta + i*dtheta)*std::cos(rubbing.tftPhi + i*dphi);
                lcDir(i)(1) = std::sin(rubbing.tftTheta + i*dtheta)*std::sin(rubbing.tftPhi + i*dphi);
                lcDir(i)(2) = std::cos(rubbing.tftTheta + i*dtheta);
            }
            epsilonr.updateEpsilonr(lcParam.epsr_para, lcParam.epsr_perp, lcDir);
        }

        ///reset LC parameters, it will change the updater
        void resetConditions(const LCParamters _lcParam);
        ///reset rubbing conditions, it will reset LC directors
        void resetConditions(const RubbingCondition _rubbing);
        ///use vector form 
        void createVectorFormUpdater(const Potential& _pot, double dt);
    private:
        size_t layerNum;
        double cellgap;
        ///LC parameters
        LCParamters lcParam;
        ///rubbing conditions
        RubbingCondition rubbing;
        ///LC directors
        DIRVEC lcDir;
        ///epsilon values
        Epsilon& epsilonr;
        ///The functor used to update LC directors
        std::shared_ptr<LCSolverBase> lcUpdater;
    };

    ///only for potentials in LC, potentials between pi and LC will be calculated with capacitance.
    class Potential{
    public:
        Potential(size_t _size, const Epsilon& _epsilonr, double _dz):epsilonr(_epsilonr), dz(_dz){
            potential.resize(_size);
            for(auto& i: potential) i = 0.0;
            EFieldForLC.resize(_size);
            for(auto& i: EFieldForLC) i = 0.0;
        }
        const size_t getSize(){return potential.size();}
        DOUBLEARRAY1D getPotentials()const{return potential;}
        const double& getPotentials(size_t i)const{return potential[i];}
        const DOUBLEARRAY1D& getEfieldForLC()const{return EFieldForLC;}
        ///enter the voltage for static calculation and enter time (int ms) for dynamic calculation.
        void update(double);
        void createDynamicUpdatePolicy(std::shared_ptr<LCD::WaveformBase> _voltWavePtr = std::shared_ptr<LCD::WaveformBase>());
        void createStaticUpdatePolicy();
    private:
        const double dz;
        ///epsilons
        const Epsilon& epsilonr;
        /// potnetials inside LC
        DOUBLEARRAY1D potential;
        /// electric fields on grid points for LC update, arranged in the same index with LC directors.
        DOUBLEARRAY1D EFieldForLC;
        /// The functor used to update potentials
        std::shared_ptr<PotentialCalculate> potentialUpdater;
    };


    class PotentialCalculate{
    public:
        PotentialCalculate(DOUBLEARRAY1D& _pot, DOUBLEARRAY1D& _EFieldForLC, const Epsilon& _epsilons, double _dz);
        virtual void update(double) = 0;
    protected:
        ///give voltage B.C. on the LC layer
        void calculate(double volt);
        double dz{0.0};
        DOUBLEARRAY1D& pot;
        DOUBLEARRAY1D& EFieldForLC;
        const Epsilon& epsilonr;
        ///This is for solving the matrix. In AX= b, it's the A.
        Eigen::MatrixX3d matrixX3d;
        ///This is for solving the matrix. In AX= b, it's the b.
        Eigen::VectorXd b;
        DOUBLEARRAY1D x;
    };

    /**
    This class serves as UpdatePolicy of Potential to update potentials and B.C. can be changed through a member function.
    */
    class PotentialSolversForStatic: public PotentialCalculate{
    public:
        PotentialSolversForStatic(DOUBLEARRAY1D& _pot, DOUBLEARRAY1D& _EFieldForLC, const Epsilon& _epsilons, double _dz);
        ///input voltage on the cell boundary, it will be transformed to the voltage BC. on the LC layer.
        virtual void update(double volt);
    };

    /**
    This class serves as UpdatePolicy of Potential to update potentials and it uses a Waveform objects to 
    decides BC.
    */
    class PotentialSolversForDynamic: public PotentialCalculate{
    public:
        PotentialSolversForDynamic(DOUBLEARRAY1D& _pot, DOUBLEARRAY1D& _EFieldForLC, const Epsilon& _epsilons, double _dz,
            std::shared_ptr<LCD::WaveformBase> _voltWavePtr = std::shared_ptr<LCD::WaveformBase>());
        virtual void update(double t);
    protected:
        std::shared_ptr<LCD::WaveformBase> voltWavePtr;
    };

    class LCSolverBase{
    public:
        void changeLCParams(LCParamters _param){
            k11 = _param.k11;
            k22 = _param.k22;
            k33 = _param.k33;
            epsr_para = _param.epsr_para;
            epsr_perp = _param.epsr_perp;
            delta_epsr = epsr_para - epsr_perp;
            gamma = _param.gamma;
            q0 = _param.q0;
        }
        void change_dz(double _dz){dz = _dz;}
        void change_dt(double _dt){dt = _dt;}
        virtual double update() = 0;
    protected:
        LCSolverBase(const Potential& _pot, LCDirector& _lcDir, Epsilon& _epsilonr, LCParamters _lcParam
            , double _dz, double _dt);
        void updateEpsilonr();
        double dz;
        double dt;
        double k11;
        double k22;
        double k33;
        double epsr_para;
        double epsr_perp;
        double delta_epsr;
        double gamma;
        double q0;
        /// tempararily put new directors.
        DIRVEC tempDir;
        ///LC parameters
        LCParamters lcParam;
        const Potential& pot;
        LCDirector& lcDir;
        Epsilon& epsilonr;
    };

    /**
    This class serves as UpdatePolicy of LCDirector to update LC directors one step forward with the vetor form calculation.
    */
    class LCVecUpdater:public LCSolverBase{
    public:
        LCVecUpdater(const Potential& _pot, LCDirector& _lcDir, Epsilon& _epsilonr, LCParamters _lcParam, double _dz, double _dt);
        ///return the residual, it is defined the maximum error of one director component among directors
        virtual double update();
    protected:
        /// temp data for variation calculations.
        DIRVEC temp;
    };
};

#endif
