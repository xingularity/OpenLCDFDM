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

#include "LCD_Optics2x2.hpp"

using namespace LCDOptics;

Optical2X2OneLayer<ISOType>::Optical2X2OneLayer(double thickness, NKData _nk, OpticalMaterialClass _layerMaterialClass):
Optical2X2OneLayerBase(_layerMaterialClass), nk_in(_nk){
    materialtype = ISOType;
    this->d = thickness;
    nk = _nk;
}

void Optical2X2OneLayer<ISOType>::interpolateNKForLambdas(DOUBLEARRAY1D lambdas_){
    nkInterpolator.resetTargetLambdas(lambdas_);
    nkInterpolator.interpolate(nk_in, nk);
}

COMPD Optical2X2OneLayer<ISOType>::findNK(double lambda){
    NKData::iterator iter = nk.find(lambda);
    //if find in map, then return. Otherwise, interposlate immediately
    return (iter != nk.end())?iter->second:nkInterpolator.interpolate(nk, lambda);
}

double Optical2X2OneLayer<ISOType>::calcJonesMatrix(JONESMAT& M, Angle& iang, double lambda, double lastn){
    double ts=0.0, tp=0.0;
    double nAvg = findNK(lambda).real();
    double theta_i, phi;
    std::tie(theta_i, phi) = iang;
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
    double theta_i, phi;
    std::tie(theta_i, phi) = iang;
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

Optical2X2OneLayer<UniaxialType>::Optical2X2OneLayer(double thickness, NKoNKeData _nk, OpticalMaterialClass _layerMaterialClass):
Optical2X2OneLayerBase(_layerMaterialClass), nk_in(_nk){
    materialtype = UniaxialType;
    this->d = thickness;
    nk = _nk;
}

void Optical2X2OneLayer<UniaxialType>::interpolateNKForLambdas(DOUBLEARRAY1D lambdas_){
    nkInterpolator.resetTargetLambdas(lambdas_);
    nkInterpolator.interpolate(nk_in, nk);
}

NKoNKe Optical2X2OneLayer<UniaxialType>::findNK(double lambda){
    NKoNKeData::iterator iter = nk.find(lambda);
    //if find in map, then return. Otherwise, interpolate immediately
    return (iter != nk.end())?iter->second:nkInterpolator.interpolate(nk, lambda);
}

double Optical2X2OneLayer<UniaxialType>::calcJonesMatrix(JONESMAT& M, Angle& iang, double lambda, double lastn){
    if (axisVec.size() < 1) throw runtime_error("No optical axis data! reported in Optical2X2OneLayer<UniaxialType>::calcJonesMatrix");
    NKoNKe nn = findNK(lambda);
    COMPD no, ne;
    std::tie(no, ne)=nn;
    //double nAvg = no.real()*2.0/3.0 + ne.real()/3.0;
    double nAvg = no.real();
    double theta_i, phi;
    std::tie(theta_i, phi) = iang;

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
    COMPD no, ne;
    std::tie(no, ne)=nn;
    //double nAvg = no.real()*2.0/3.0 + ne.real()/3.0;
    double nAvg = no.real();

    double theta_i, phi;
    std::tie(theta_i, phi) = iang;
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
