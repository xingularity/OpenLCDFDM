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

#include "LCD1D_ExtendedJones.hpp"
#include <stdexcept>
#include <algorithm>
#include <typeinfo>
#include <omp.h>

using namespace LCD1D;

ExtendedJones::ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles,
    const double targetLambda, LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifLambertian, bool _ifStokes):
ExtendedJonesBase(_materials, _inAngles, targetLambda, lightSrcSpectrum_)
{
    ifLambertian = _ifLambertian;
    findLCLayerInMaterialList();
    if (_ifStokes)checkIfCalcStokes();
    if (ifCalcStokes)resetStokes();
}

ExtendedJones::ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
const double end_lambda_, const double step_lambda_, LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifLambertian, bool _ifStokes):
ExtendedJonesBase(_materials, _inAngles, start_lambda_, end_lambda_, step_lambda_, lightSrcSpectrum_)
{
    ifLambertian = _ifLambertian;
    findLCLayerInMaterialList();
    if (_ifStokes)checkIfCalcStokes();
    if (ifCalcStokes)resetStokes();
}

void ExtendedJones::checkIfCalcStokes(){
    ifCalcStokes = false;
    //There must be at least one polarizer layer.
    //check number of polarizers first
    if (polarizerLayersIndex.size() < 1){
        std::cout << "Number of the polarizer layer < 1, no Stokes calculation" << std::endl;
        return;
    }
    /*
    int index = polarizerLayersIndex[0];
    std::shared_ptr<Optical2X2OneLayer<UniaxialType> > tempLayerPtr;
    tempLayerPtr = std::dynamic_pointer_cast<Optical2X2OneLayer<UniaxialType>,Optical2X2OneLayerBase>(matLayers[index]);
    if (!tempLayerPtr) throw runtime_error("error type conversion while trying to find optical axis of the first polarizers.");
    DIRVEC axes = tempLayerPtr -> getAxes();
    //polarizer's optical axis should have almost 90 degree angle
    BzVECD3D axis = axes(0);
    if ((std::abs(acos(axis(2))) - M_PI/2.0) > 1.0e-10) && ((std::abs(acos(axis(2))) - 3.0*M_PI/2.0) > 1.0e-10){
        std::cout << "The axis of the first polarizer doesn't approach to 90 or 270 angle." << std::endl;
        return;
    }
    */
    std::cout << "Stokes calculation will be used." <<std::endl;
    ifCalcStokes = true;
    return;
}

void ExtendedJones::calculateExtendedJones(){
    //reset before calculation to be sure.
    resetTransmissions();
    resetTransEachLambda();
    if (lambdas.size() == 0)
        throw std::runtime_error("can't calculate extended Jones matrix without any given wavelength");
    //interpolate all NK data to target wavelengths
    for (int i = 0; i < matLayers.size(); ++i)
        matLayers[i]->interpolateNKForLambdas(lambdas);
    if (ifCalcStokes){
        resetStokes();
        throw std::runtime_error("Not yet support Stokes in ExtendedJones::calculateExtendedJones()");
    }else{
        if (lambdas.size() == 1){
            calculateOneLambdaNoStokes(0);
            //Lambertian light source
            if (ifLambertian){
                for(int i = 0; i < transmissions.size(); ++i)
                    for(int j = 0; j < transmissions[i].size(); ++j)
                        transEachLambda[0][i][j]*=cos(std::get<0>(inAngles[i][j]));
            }
            transmissions = transEachLambda[0];
        }
        else{
            double denominator = 0.0;
            double step_lambda = lambdas[1] - lambdas[0];
            for (int i =0; i < lightSourceSpectrum.size(); ++i)
                denominator += lightSourceSpectrum[i] * yBarOfLambda[i] * step_lambda;

            #pragma omp parallel for
            for (int i = 0; i < lambdas.size();++i){
                std::cout << omp_get_num_threads() <<std::endl;
                calculateOneLambdaNoStokes(i);
            }

            //Lambertian light source
            if (ifLambertian){
                for (int i = 0; i < lambdas.size();++i)
                    for(int j = 0; j < transmissions.size(); ++j)
                        for(int k = 0; k < transmissions[i].size(); ++k)
                            transEachLambda[i][j][k]*=cos(std::get<0>(inAngles[j][k]));
            }

            //calculate multi-wavelength transmissions
            for (int i = 0; i < lambdas.size();++i)
                for(int j = 0; j < transmissions.size(); ++j)
                    for(int k = 0; k < transmissions[j].size(); ++k)
                        transmissions[j][k] += transEachLambda[i][j][k] * lightSourceSpectrum[i] * yBarOfLambda[i] * step_lambda;

            for(int i = 0; i < transmissions.size(); ++i)
                for(int j = 0; j < transmissions[i].size(); ++j)
                    transmissions[i][j]/=denominator;
        }
    }
}

void ExtendedJones::calculateOneLambdaNoStokes(int iLambda){
    //assuming incident from air
    double lastn = nAir;
    double ts = 1.0, tp = 1.0;
    double lambda = lambdas[iLambda];
    for (int i = 0; i < inAngles.size(); ++i)
        for(int j = 0; j < inAngles[i].size(); ++j){
            JONESMAT M;
            M << 1.0,0.0,0.0,1.0;
            Angle inAngle = inAngles[i][j];
            lastn = nAir;
            for (int k =0; k < matLayerNum; k++){
                lastn = matLayers[k]->calcJonesMatrix(M, inAngle, lambda, lastn);
            }
                const double& theta_i = std::get<0>(inAngle);
                  //refraction back to air
                const double& theta_r = std::get<0>(inAngles[i][j]);
                LCDOptics::EigenM22 tMat;
                tMat << (2.0*lastn*cos(theta_i)/(lastn*cos(theta_i)+nAir*cos(theta_r))),0,0,
                        (2.0*lastn*cos(theta_i)/(lastn*cos(theta_r)+nAir*cos(theta_i)));
                M=tMat*M;
                  //put tramissions into transEachLambda first
                transEachLambda[iLambda][i][j]=0.5*(pow(std::abs(M(0,0)),2.0)+pow(std::abs(M(0,1)),2.0)
                +pow(std::abs(M(1,0)),2.0)+pow(std::abs(M(1,1)),2.0));
        }
}

void ExtendedJones::calculateOneLambdaWithStokes(int iLambda){
    //If the program comes here, it means there is at least one polarizer and its theta angle approaches to 90 or 270 degree.
    //assuming incident from air
    resetTransEachLambda();
    double lastn = nAir;
    double ts = 1.0, tp = 1.0;
    double lambda = lambdas[iLambda];
    unsigned int polarCalcStart = *(std::min_element(polarizerLayersIndex.begin(), polarizerLayersIndex.end()));
    //polarCalcEnd points to the next layer of the final polarizer.
    unsigned int polarCalcEnd = (polarizerLayersIndex.size() > 1) ? *(std::max_element(polarizerLayersIndex.begin(), polarizerLayersIndex.end())) + 1: matLayerNum;
    for (int i = 0; i < inAngles.size(); ++i)
        for(int j = 0; j < inAngles[i].size(); ++j){
            JONESMAT M;
            M << 1.0,0.0,0.0,1.0;
            Angle inAngle = inAngles[i][j];
            POLARTRACE lightPolar(0); //no elements at init
            lastn = nAir;
            for (int k =0; k < matLayerNum; k++){
                if ((k >= polarCalcStart) && (k < polarCalcEnd))
                    lastn = matLayers[k]->calcJonesMatrix(M, lightPolar, inAngle, lambda, lastn);
                else
                    lastn = matLayers[k]->calcJonesMatrix(M, inAngle, lambda, lastn);
            }
            const double& theta_i = std::get<0>(inAngle);
              //refraction back to air
            const double& theta_r = std::get<0>(inAngles[i][j]);
            LCDOptics::EigenM22 tMat;
            tMat << (2.0*lastn*cos(theta_i)/(lastn*cos(theta_i)+nAir*cos(theta_r))),0,0,
                    (2.0*lastn*cos(theta_i)/(lastn*cos(theta_r)+nAir*cos(theta_i)));
            M=tMat*M;
              //put tramissions into transEachLambda first
            transEachLambda[iLambda][i][j]=0.5*(pow(std::abs(M(0,0)),2.0)+pow(std::abs(M(0,1)),2.0)
            +pow(std::abs(M(1,0)),2.0)+pow(std::abs(M(1,1)),2.0));

            Eigen::Vector2cd polar = lightPolar.back();
            polar = tMat*polar;
            lightPolar.push_back(polar);

            //start to convert polarization to Stokes vector with assumtion that s0 = 1
            STOKESTRACE stokesTrace;
            for (unsigned int k =0; k < lightPolar.size(); ++k){
                Eigen::Vector2cd polar = lightPolar[k];
                Eigen::Vector3d stokeVector;
                double psi = std::atan2(polar[1].real(), polar[0].real());
                double delta = std::arg(polar[1] / polar[0]);
                stokeVector << std::cos(2.0*psi), std::sin(2.0*psi)*std::cos(delta), std::sin(2.0*psi)*std::sin(delta);
                stokesTrace.push_back(stokeVector);
            }
            stokes[iLambda][i][j] = stokesTrace;
        }
}

const TRANSRESULT& ExtendedJones::getTransmissions(){
    return transmissions;
}

const STOKESRESULT& ExtendedJones::getStokes(){
    return stokes;
}

void ExtendedJones::resetTransmissions(){
    transmissions.clear();
    transmissions = DOUBLEARRAY2D(inAngles.size(), DOUBLEARRAY1D(inAngles[0].size(), 0.0));
}

void ExtendedJones::resetTransEachLambda(){
    transEachLambda.clear();
    transEachLambda = std::vector<TRANSRESULT>(lambdas.size(), DOUBLEARRAY2D(inAngles.size(), DOUBLEARRAY1D(inAngles[0].size(), 0.0)));
}
void ExtendedJones::resetStokes(){
    stokes.clear();
    stokes = STOKESRESULT(lambdas.size(), std::vector<std::vector<STOKESTRACE> >(inAngles.size(),
        std::vector<STOKESTRACE>(inAngles[0].size())));
}

void ExtendedJones::resetLCDiretors(DIRVEC _in){
    if (lcLayerindex >= 0){
        LCDOptics::Optical2x2UnixialPtr tempLayerPtr;
        tempLayerPtr = std::dynamic_pointer_cast<LCDOptics::Optical2X2OneLayer<LCDOptics::UniaxialType>,
            LCDOptics::Optical2X2OneLayerBase>(matLayers[lcLayerindex]);
        if (tempLayerPtr) tempLayerPtr->resetAxes(_in);
    }
}
