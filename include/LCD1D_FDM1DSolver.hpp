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
        Epsilon(size_t _lc_layernum, double _dz){lcLayerNum=_lc_layernum; epsr33.resize(lcLayerNum); dz=_dz}
        ///set TFTPI parameters
        void setTFTPI(DielecParameters _param){tftpi_epsr = _param.epsr; tftpiThick = _param.thick;}
        ///set CFPI parameters
        void setCFPI(DielecParameters _param){cfpi_epsr = _param.epsr; cfpiThick = _param.thick;}
        ///get epsilon distribution in LC
        DOUBLEARRAY1D getLCEpsr33()const{return epsr33;}
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
        void calculateCapacitance() const;
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
    template <class UpdatePolicy>
    class LCDirector{
        ///C++11 syntax allowed
        friend UpdatePolicy;
    public:
        ///Constructor
        LCDirector(size_t _lcLayerNum, const LCParamters _lcParam, const RubbingCondition _rubbing, 
            Epsilon& _epsilionr):lcParam(_lcParam), rubbing(_rubbing), epsilonr(_epsilonr){
            
            lcDir.resizeAndPreserve(_lcLayerNum+1);
            this->resetDirectors();
        }
        size_t getSize()const{return dir.extent(0);}
        size_t getLayerNum()const{return dir.extent(0) - 1;}
        const DIRVEC getDirectors()const {return lcDir;}
        LCParamters getLCParams()const {return lcParam;}
        RubbingCondition getLCRubbing()const {return rubbing;}
        double getDirectors(int iz, int i)const {return lcDir(iz)(i);}
        const BzVECD3D getDirectors(int iz)const {return lcDir(iz);}
        ///reset LC directors according to the rubbing BC.
        void resetLCDirectors(){
            size_t layerNum = lcDir.extent(0)-1;
            double dtheta = (lcParam.cfTheta - lcParam.tftTheta)/layerNum;
            double dphi = lcParam.totalTwis/layerNum;
            for (int i = 0; i < layerNum + 1; ++i){
                lcDir(i)(0) = std::sin(rubbing.tftTheta + i*dtheta)*std::cos(rubbing.tftPhi + i*dphi);
                lcDir(i)(1) = std::sin(rubbing.tftTheta + i*dtheta)*std::sin(rubbing.tftPhi + i*dphi);
                lcDir(i)(2) = std::cos(rubbing.tftTheta + i*dtheta);
            }
            epeilonr.updateEpsilonr(lcParam.epsr_para, lcParam.epsr_perp, lcDir);
        }
        ///reset LC parameters
        void resetConditions(const LCParamters _lcParam){lcParam = _lcParam;}
        ///reset rubbing conditions
        void resetConditions(const RubbingCondition _rubbing){rubbing = _rubbing;}
    private:
        ///LC parameters
        LCParamters lcParam;
        ///rubbing conditions
        RubbingCondition rubbing;
        ///LC directors
        DIRVEC lcDir;
        ///epsilon values
        Epsilon& epsilonr;
        ///The functor used to update LC directors
        std::shared_ptr<UpdatePolicy> lcUpdater;
    };

    ///only for potentials in LC, potentials between pi and LC will be calculated with capacitance.
    template <class UpdatePolicy>
    class Potential{
        friend UpdatePolicy;
    public:
        Potential(size_t _size, const Epsilon& _epsilonr):epsilonr(_epsilonr){
            potential.resize(_size);
            for(auto& i: potential) i = 0.0;
        }
        const getSize(){return potential.size();}
        DOUBLEARRAY1D getPotentials()const{return potential;}
        const double& getPotentials(size_t i)const{return potential[i];}
        const DOUBLEARRAY1D& getEfieldForLC()const{return EFieldForLC;}
        ///enter the voltage for static calculation and enter time (int ms) for dynamic calculation.
        void update(double);
    private:
        ///epsilons
        const Epsilon& epsilonr;
        /// potnetials inside LC
        DOUBLEARRAY1D potential;
        /// electric fields on grid points for LC update, arranged in the same index with LC directors.
        DOUBLEARRAY1D EFieldForLC;
        /// The functor used to update potentials
        std::shared_ptr<UpdatePolicy> potentialUpdater;
    };
    
    template<class UpdatePolicy>
    Potential<UpdatePolicy>::update(double anything){
        throw runtime_error("This function can't be called, one must use specialized version.");
    }

    template<>
    Potential<PotentialSolversForStatic>::update(double volt){
        potentialUpdater->update(volt);
    }

    template<>
    Potential<PotentialSolversForDynamic>::update(double t){
        potentialUpdater->update(t);
    }


    /**
    This class serves as UpdatePolicy of Potential to update potentials and B.C. can be changed through a member function.
    */
    class PotentialSolversForStatic{
    public:
        PotentialSolver1D(Potnetial& _pot, const Epsilon& _epsilons, double _dz);
        ///input voltage on the cell boundary, it will be transformed to the voltage BC. on the LC layer.
        void update(double volt);
    protected:
        double dz{0.0};
        Potential& potentials;
        const Epsilon& epsilonr;
    };

    /**
    This class serves as UpdatePolicy of Potential to update potentials and it uses a Waveform objects to 
    decides BC.
    */
    class PotentialSolversForDynamic{
    public:
        PotentialSolver1D(Potnetial& _pot, LCDirector& _lcDir, double _dz, DielecParameters _tftpi,
            std::shared_ptr<LCD::WaveformBase> _voltWavePtr = std::shared_ptr<LCD::WaveformBase>());
        void update(double t);
    protected:
        double dz{0.0};
        Potential& potentials;
        Epsilon& epsilonr;
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
    protected:
        LCSolverBase(const Potnetial& _pot, LCDirector& _lcDir, Epsilon& _epsilonr, LCParamters _lcParam
            , double _dz, double _dt):epsilonr(_epsilonr), lcParam(_lcParam), potentials(_pot), lcDir(_lcDir), dz(_dz), dt(_dt){
                tempDir.resizeAndPreserve(lcDir.getSize());
            }
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
        const Potential& potentials;
        LCDirector& lcDir;
        Epsilon& epsilonr;
    };

    /**
    This class serves as UpdatePolicy of LCDirector to update LC directors one step forward with the vetor form calculation.
    */
    class LCVecUpdater:public LCSolverBase{
    public:
        LCSolver(const Potnetial& _pot, LCDirector& _lcDir, LCParamters _lcParam, double _dz, double _dt):
            LCSolverBase(_pot, _lcDir, _lcParam, _dz, _dt){
                temp.resizeAndPreserve(lcDir.getSize());
            }
        ///return the residual, it is defined the maximum error of one director component among directors
        double update();
    };
    protected:
        /// temp data for variation calculations.
        DIRVEC temp;
};

#endif
