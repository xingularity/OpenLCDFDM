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

#ifndef LCD_OPTICS2X2_HPP
#define LCD_OPTICS2X2_HPP

#include "LCD_OpticsBaseDef.hpp"

namespace LCDOptics{
    using LCD::EigenC22;
    using LCD::EigenM22;
    using Eigen::Vector3d;
    using POLARTRACE = std::vector<Eigen::Vector2cd>;
    ///Jones matrix is a 2x2 matrix.
    typedef EigenC22 JONESMAT;
    enum OpticalMaterialClass{
        OPT_GLASS = 0,
        OPT_ISOTROPIC = 1,
        OPT_UNIAXIAL = 2,
        OPT_LCMATERIAL = 3,
        OPT_POLARIZER = 4,
        OPT_BIAXIAL = 5
    };

    /**
    Optical2X2OneLayerBase represents base class of one layer of material when calculating 2X2 extended Jones matrix.<br/>
    */
    class Optical2X2OneLayerBase{
    public:
        Optical2X2OneLayerBase(OpticalMaterialClass _layerMaterialClass){layerMaterialClass = _layerMaterialClass;}
        ///calculate one jones matrix and return its average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, Angle& iang, double lambda, double lastn) = 0;
        ///calculate one jones matrix and light polarization, return the average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, POLARTRACE& lightPolar, Angle& iang, double lambda, double lastn) = 0;
        virtual void interpolateNKForLambdas(DOUBLEARRAY1D lambdas_) = 0;
        OpticalMaterialClass opticalLayerKind(){return layerMaterialClass;}
        const DIRVEC getAxes(){return axisVec;}
    protected:
        ///calculate polarizations based on input Jones matrix
        void calculatePolarization(const EigenC22& m, POLARTRACE& lightPolar){
            if (lightPolar.size() == 0) return;
            Eigen::Vector2cd polar = lightPolar.back();
            polar = m*polar;
            lightPolar.push_back(polar);
        }
        ///calculate polarizations based on input Jones matrix
        void calculatePolarization(const EigenM22& m, POLARTRACE& lightPolar){
            if (lightPolar.size() == 0) return;
            Eigen::Vector2cd polar = lightPolar.back();
            polar = m*polar;
            lightPolar.push_back(polar);
        }
        ///thickness of this layer
        double d;
        ///material type of this layer
        OpticMaterialType materialtype;
        ///How may sublayer in the layer.
        size_t layernum;
        ///data of optical axis
        DIRVEC axisVec;
        OpticalMaterialClass layerMaterialClass;
    };

    template<int E>
    class Optical2X2OneLayer: public std::false_type{};

    /**
    Optical2X2IsotropicLayer represents base class of one layer of isotropic material when calculating 2X2 extended Jones matrix.<br/>
    <em>It calculates incident tstp matrix, but doesn't do output tstp matrix.</em><br/>
    */
    template<>
    class Optical2X2OneLayer<ISOType>: public Optical2X2OneLayerBase{
    public:
        ///For isotropic material.
        Optical2X2OneLayer(double thickness, NKData _nk, OpticalMaterialClass _layerMaterialClass = OPT_ISOTROPIC);
        ///calculate one jones matrix and return its average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, Angle& iang, double lambda, double lastn);
        ///calculate one jones matrix and light polarization, return the average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, POLARTRACE& lightPolar, Angle& iang, double lambda, double lastn);
        virtual void interpolateNKForLambdas(DOUBLEARRAY1D lambdas_);
    private:
        ///find nk for corresponding lambda
        COMPD findNK(double lambda);
        ///nk dispersion for isotropic material
        NKData nk;
        ///Interpolation of nk distribution
        SpectrumInterpolator<NKData> nkInterpolator;
    };

    Optical2X2OneLayer<ISOType>::Optical2X2OneLayer(double thickness, NKData _nk, OpticalMaterialClass _layerMaterialClass):
    Optical2X2OneLayerBase(_layerMaterialClass){
        materialtype = ISOType;
        this->d = thickness;
        nk = _nk;
    }

    void Optical2X2OneLayer<ISOType>::interpolateNKForLambdas(DOUBLEARRAY1D lambdas_){
        nkInterpolator.resetTargetLambdas(lambdas_);
        NKData new_nk;
        nkInterpolator.interpolate(nk, new_nk);
        nk = new_nk;
    }

    COMPD Optical2X2OneLayer<ISOType>::findNK(double lambda){
        NKData::iterator iter = nk.find(lambda);
        //if find in map, then return. Otherwise, interposlate immediately
        return (iter != nk.end())?iter->second:nkInterpolator.interpolate(nk, lambda);
    }

    double Optical2X2OneLayer<ISOType>::calcJonesMatrix(JONESMAT& M, Angle& iang, double lambda, double lastn){
        double ts=0.0, tp=0.0;
        double nAvg = findNK(lambda).real();
        const double& theta_i = std::get<0>(iang);
        const double& phi = std::get<1>(iang);
        //refraction angle
        double theta_r = std::asin(lastn*sin(theta_i)/nAvg);
        EigenM22 tMat;
        tMat << 0,0,0,0;
        ts=2.0*lastn*cos(theta_i)/(lastn*cos(theta_i)+nAvg*cos(theta_r));
        tp=2.0*lastn*cos(theta_i)/(lastn*cos(theta_r)+nAvg*cos(theta_i));
        tMat(0,0)=ts;tMat(1,1)=tp;
        M=tMat*M;
        std::get<0>(iang)=theta_r;
        return nAvg;
    }

    double Optical2X2OneLayer<ISOType>::calcJonesMatrix(JONESMAT& M, POLARTRACE& lightPolar,
    Angle& iang, double lambda, double lastn){
        double ts=0.0, tp=0.0;
        double nAvg = findNK(lambda).real();
        const double& theta_i = std::get<0>(iang);
        const double& phi = std::get<1>(iang);
        //refraction angle
        double theta_r = std::asin(lastn*sin(theta_i)/nAvg);
        EigenM22 tMat;
        tMat << 0,0,0,0;
        ts=2.0*lastn*cos(theta_i)/(lastn*cos(theta_i)+nAvg*cos(theta_r));
        tp=2.0*lastn*cos(theta_i)/(lastn*cos(theta_r)+nAvg*cos(theta_i));
        tMat(0,0)=ts;tMat(1,1)=tp;
        M=tMat*M;
        //calculate light polarization
        calculatePolarization(tMat, lightPolar);

        std::get<0>(iang)=theta_r;
        return nAvg;
    }

    /**
    Optical2X2UniaxialLayer represents base class of one layer of uniaxial material when calculating 2X2 extended Jones matrix.<br/>
    <em>It calculates incident tstp matrix, but doesn't do output tstp matrix.</em><br/>
    <em>Notice that One must reset optical axis before every calculation</em><br/>
    */
    template<>
    class Optical2X2OneLayer<UniaxialType>: public Optical2X2OneLayerBase{
    public:
        ///For uniaxial material.
        Optical2X2OneLayer(double thickness, NKoNKeData _nk, OpticalMaterialClass _layerMaterialClass = OPT_UNIAXIAL);
        void resetDirectors(DIRVEC _in);
        ///calculate one jones matrix and return its average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, Angle& iang, double lambda, double lastn);
        ///calculate one jones matrix and light polarization, return the average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, POLARTRACE& lightPolar, Angle& iang, double lambda, double lastn);
        virtual void interpolateNKForLambdas(DOUBLEARRAY1D lambdas_);
    private:
        ///find nk for corresponding lambda. If find in map, then return. Otherwise, interpolate immediately
        NKoNKe findNK(double lambda);
        ///nk dispersion for uniaxial material
        NKoNKeData nk;
        ///Interpolation of nk distribution
        SpectrumInterpolator<NKoNKeData> nkInterpolator;
    };

    Optical2X2OneLayer<UniaxialType>::Optical2X2OneLayer(double thickness, NKoNKeData _nk, OpticalMaterialClass _layerMaterialClass):
    Optical2X2OneLayerBase(_layerMaterialClass){
        materialtype = UniaxialType;
        this->d = thickness;
        nk = _nk;
    }

    void Optical2X2OneLayer<UniaxialType>::interpolateNKForLambdas(DOUBLEARRAY1D lambdas_){
        nkInterpolator.resetTargetLambdas(lambdas_);
        NKoNKeData new_nk;
        nkInterpolator.interpolate(nk, new_nk);
        nk = new_nk;
    }

    NKoNKe Optical2X2OneLayer<UniaxialType>::findNK(double lambda){
        NKoNKeData::iterator iter = nk.find(lambda);
        //if find in map, then return. Otherwise, interpolate immediately
        return (iter != nk.end())?iter->second:nkInterpolator.interpolate(nk, lambda);
    }

    double Optical2X2OneLayer<UniaxialType>::calcJonesMatrix(JONESMAT& M, Angle& iang, double lambda, double lastn){
        NKoNKe nn = findNK(lambda);
        const COMPD& no = std::get<0>(nn);
        const COMPD& ne = std::get<1>(nn);
        double nAvg = no.real()*2.0/3.0 + ne.real()/3.0;

        const double& theta_i = std::get<0>(iang);
        const double& phi = std::get<1>(iang);
        double theta_r=asin(lastn*sin(theta_i)/nAvg);

        //tangential part of k, reference to Dr. Pochi Yeh's book
        double alpha=(2.0*M_PI/lambda)*lastn*sin(theta_i)*cos(phi);
        double beta=(2.0*M_PI/lambda)*lastn*sin(theta_i)*sin(phi);

        EigenM22 tMat; tMat << 0,0,0,0;
        double ts=0.0, tp=0.0;
        ts=2.0*lastn*cos(theta_i)/(lastn*cos(theta_i)+nAvg*cos(theta_r));
        tp=2.0*lastn*cos(theta_i)/(lastn*cos(theta_r)+nAvg*cos(theta_i));
        tMat(0,0)=ts;tMat(1,1)=tp;

        //Multiply ts tp matrix. M is now going into the zero thickness imaginary layer
        M=tMat*M;

        Vector3d zz; zz <<0.0,0.0,1.0;
        Vector3d kr; kr << sin(theta_r)*cos(phi), sin(theta_r)*sin(phi), cos(theta_r);
        Vector3d s,p;
        if (kr==zz){
            s << 0.0,-1.0,0.0;
            p << 1.0,0.0,0.0;
        }
        else{
            s=kr.cross(zz);s.normalize();
            p=kr.cross(s);p.normalize();
        }
        Vector3d caxis_temp; caxis_temp << axisVec(0)(0), axisVec(0)(1), axisVec(0)(2);
        Vector3d o=kr.cross(caxis_temp);o.normalize();
        Vector3d e=o.cross(kr);e.normalize();
        EigenM22 sp_eo;
        sp_eo << s.dot(e), p.dot(e), s.dot(o), p.dot(o);
        //convert sp coordinate to eo coordinate.
        M=sp_eo*M;

        double ld=d/layernum;
        double angle[2]; //temp array for theta and phi of one axis
        for (int i=0; i< layernum;i++){
        //calcutate kez, koz
            angle[0] = acos(axisVec(i)(2));
            angle[1] = atan2(axisVec(i)(1), axisVec(i)(0));

            double kd=alpha*cos(angle[1])+beta*sin(angle[1]);
            double keb=alpha*sin(angle[1])*(-1.0)+beta*cos(angle[1]);
            COMPD u=pow(sin(angle[0])/ne,2.0)+pow(cos(angle[0])/no,2.0);
            COMPD v=kd*sin(2.0*angle[0])*(pow(1.0/ne,2.0)-pow(1.0/no,2.0));
            COMPD w=(pow(kd*cos(angle[0]),2.0)+pow(keb,2.0))/pow(ne,2.0)
                      +pow(kd*sin(angle[0])/no,2.0)
                      -(2.0*M_PI/lambda)*(2.0*M_PI/lambda);
            COMPD kez=(v+sqrt(v*v-4.0*u*w))/2.0/u;
            COMPD koz=sqrt(no*no*(2.0*M_PI/lambda)*(2.0*M_PI/lambda)-alpha*alpha-beta*beta);

            M(0,0)*=exp(COMPD(0.0,-1.0)*kez*ld);
            M(0,1)*=exp(COMPD(0.0,-1.0)*kez*ld);
            M(1,0)*=exp(COMPD(0.0,-1.0)*koz*ld);
            M(1,1)*=exp(COMPD(0.0,-1.0)*koz*ld);

            if (i<layernum-1){
                caxis_temp << axisVec(i+1)(0), axisVec(i+1)(1), axisVec(i+1)(2);
                Vector3d o1=kr.cross(caxis_temp);o1.normalize();
                Vector3d e1=o1.cross(kr);e1.normalize();
                EigenM22 eo_trans;
                eo_trans << e.dot(e1), o.dot(e1), e.dot(o1), o.dot(o1);
                M=eo_trans*M;
                o=o1;e=e1;
            }
        }
        //calculate output D matrix
        EigenM22 eo_sp;
        eo_sp << e.dot(s), o.dot(s), e.dot(p), o.dot(p);
        M=eo_sp*M;
        std::get<0>(iang)=theta_r;
        return nAvg;
    }

    double Optical2X2OneLayer<UniaxialType>::calcJonesMatrix(JONESMAT& M,
    POLARTRACE& lightPolar, Angle& iang, double lambda, double lastn){
        NKoNKe nn = findNK(lambda);
        const COMPD& no = std::get<0>(nn);
        const COMPD& ne = std::get<1>(nn);
        double nAvg = no.real()*2.0/3.0 + ne.real()/3.0;

        const double& theta_i = std::get<0>(iang);
        const double& phi = std::get<1>(iang);
        double theta_r=asin(lastn*sin(theta_i)/nAvg);

        //tangential part of k, reference to Dr. Pochi Yeh's book
        double alpha=(2.0*M_PI/lambda)*lastn*sin(theta_i)*cos(phi);
        double beta=(2.0*M_PI/lambda)*lastn*sin(theta_i)*sin(phi);

        EigenM22 tMat; tMat << 0,0,0,0;
        double ts=0.0, tp=0.0;
        ts=2.0*lastn*cos(theta_i)/(lastn*cos(theta_i)+nAvg*cos(theta_r));
        tp=2.0*lastn*cos(theta_i)/(lastn*cos(theta_r)+nAvg*cos(theta_i));
        tMat(0,0)=ts;tMat(1,1)=tp;

        //Multiply ts tp matrix. M is now going into the zero thickness imaginary layer
        M=tMat*M;
        //calculate light polarization
        calculatePolarization(tMat, lightPolar);

        Vector3d zz; zz <<0.0,0.0,1.0;
        Vector3d kr; kr << sin(theta_r)*cos(phi), sin(theta_r)*sin(phi), cos(theta_r);
        Vector3d s,p;
        if (kr==zz){
            s << 0.0,-1.0,0.0;
            p << 1.0,0.0,0.0;
        }
        else{
            s=kr.cross(zz);s.normalize();
            p=kr.cross(s);p.normalize();
        }
        Vector3d caxis_temp; caxis_temp << axisVec(0)(0), axisVec(0)(1), axisVec(0)(2);
        Vector3d o=kr.cross(caxis_temp);o.normalize();
        Vector3d e=o.cross(kr);e.normalize();
        EigenM22 sp_eo;
        sp_eo << s.dot(e), p.dot(e), s.dot(o), p.dot(o);
        //convert sp coordinate to eo coordinate.
        M=sp_eo*M;
        calculatePolarization(sp_eo, lightPolar);
        double ld=d/layernum;
        double angle[2]; //temp array for theta and phi of one axis
        for (int i=0; i< layernum;i++){
        //calcutate kez, koz
            angle[0] = acos(axisVec(i)(2));
            angle[1] = atan2(axisVec(i)(1), axisVec(i)(0));

            double kd=alpha*cos(angle[1])+beta*sin(angle[1]);
            double keb=alpha*sin(angle[1])*(-1.0)+beta*cos(angle[1]);
            COMPD u=pow(sin(angle[0])/ne,2.0)+pow(cos(angle[0])/no,2.0);
            COMPD v=kd*sin(2.0*angle[0])*(pow(1.0/ne,2.0)-pow(1.0/no,2.0));
            COMPD w=(pow(kd*cos(angle[0]),2.0)+pow(keb,2.0))/pow(ne,2.0)
                      +pow(kd*sin(angle[0])/no,2.0)
                      -(2.0*M_PI/lambda)*(2.0*M_PI/lambda);
            COMPD kez=(v+sqrt(v*v-4.0*u*w))/2.0/u;
            COMPD koz=sqrt(no*no*(2.0*M_PI/lambda)*(2.0*M_PI/lambda)-alpha*alpha-beta*beta);

            EigenC22 phaseMatrix;
            phaseMatrix << exp(COMPD(0.0,-1.0)*kez*ld), exp(COMPD(0.0,-1.0)*kez*ld), exp(COMPD(0.0,-1.0)*koz*ld), exp(COMPD(0.0,-1.0)*koz*ld);
            M(0,0)*=phaseMatrix(0,0);
            M(0,1)*=phaseMatrix(0,1);
            M(1,0)*=phaseMatrix(1,0);
            M(1,1)*=phaseMatrix(1,1);
            //calculate light polarization
            calculatePolarization(phaseMatrix, lightPolar);
            if (i<layernum-1){
                caxis_temp << axisVec(i+1)(0), axisVec(i+1)(1), axisVec(i+1)(2);
                Vector3d o1=kr.cross(caxis_temp);o1.normalize();
                Vector3d e1=o1.cross(kr);e1.normalize();
                EigenC22 eo_trans;
                eo_trans << e.dot(e1), o.dot(e1), e.dot(o1), o.dot(o1);
                M=eo_trans*M;
                calculatePolarization(eo_trans, lightPolar);
                o=o1;e=e1;
            }
        }
        //calculate output D matrix
        EigenM22 eo_sp;
        eo_sp << e.dot(s), o.dot(s), e.dot(p), o.dot(p);
        M=eo_sp*M;
        if ((lightPolar.size() == 0) && (layerMaterialClass==OPT_POLARIZER)){
            //This layer is a polarizer and no stokes has been calculated yet. Then, this is the first layer to calculate Stokes.
            Eigen::Vector2cd polar;
            if (ne.imag() < no.imag()){
                //Extraordinary light has been perfectly absorped.
                polar<< 0.0, 1.0;
            }
            else{
                //Ordinary light has been perfectly absorped.
                polar << 1.0, 0.0;
            }
            //get incident polariation
            polar = eo_sp*polar;
            lightPolar.push_back(polar);
        }
        calculatePolarization(eo_sp, lightPolar);
        std::get<0>(iang)=theta_r;
        return nAvg;
    }

    void Optical2X2OneLayer<UniaxialType>::resetDirectors(DIRVEC _in){
        layernum = _in.size();
        axisVec = _in;
    }
};

#endif
