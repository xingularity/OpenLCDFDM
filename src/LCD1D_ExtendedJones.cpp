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
#include <typeinfo>

using namespace LCD1D;

ExtendedJones::ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles,const double targetLambda, LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifStokes):
ExtendedJonesBase(_materials, _inAngles, targetLambda, lightSrcSpectrum_)
{
    findLCLayerInMaterialList();
    if (_ifStokes)checkIfCalcStokes();
    if (ifCalcStokes)resetStokes();
}

ExtendedJones::ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
const double end_lambda_, const double step_lambda_, LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifStokes):
ExtendedJonesBase(_materials, _inAngles, start_lambda_, end_lambda_, step_lambda_, lightSrcSpectrum_)
{
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
    std::cout << "Stokes calculation will be used." <<std::endl;
    ifCalcStokes = true;
    return;
}

void ExtendedJones::calculateExtendedJones(){
    //reset before calculation to be sure.
    resetTransmissions();
    resetTransTemp();
    if (lambdas.size() == 0)
        throw runtime_error("can't calculate extended jones matrix without any given wavelength");
    //interpolate all NK data to target wavelengths
    for (int i = 0; i < matLayers.size(); ++i)
        matLayers[lcLayerindex]->interpolateNKForLambdas(lambdas);
    if (ifCalcStokes){

    }else{
        if (lambdas.size() == 1){
            calculateOneLambdaNoStokes(0);
            transmissions = transTemp;
        }
        else{
            double denominator = 0.0;
            double step_lambda = lambdas[1] - lambdas[0];
            for (int i =0; i < lightSourceSpectrum.size(); ++i)
                denominator += lightSourceSpectrum[i] * spectralEfficiency[i] * step_lambda;
            for (int i = 0; i < lambdas.size()){
                calculateOneLambdaNoStokes(i);
                for(int j = 0; j < transmissions.size(); ++j)
                    for(int k = 0; k < transmissions.size(); ++k)
                        transmissions[j][k] += transTemp[j][k] * lightSourceSpectrum[i] * spectralEfficiency[i] * step_lambda;

            }
            for(int i = 0; i < transmissions.size(); ++i)
                for(int j = 0; j < transmissions.size(); ++j)
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
            for (int k =0; k < matLayerNum; k++)
                lastn = matLayers.calcJonesMatrix(M, inAngle, lambda, lastn);
                const double& theta_i = std::get<0>(inAngle);
                  //refraction back to air
                const double& theta_r = std::get<0>(inAngles[i][j]);
                EigenM22 tMat;
                tMat << (2.0*lastn*cos(theta_i)/(lastn*cos(theta_i)+nAvg*cos(theta_r))),0,0,
                        (2.0*lastn*cos(theta_i)/(lastn*cos(theta_r)+nAvg*cos(theta_i)));
                M=tMat*M;
                  //put tramissions into transTemp first
                transTemp[i][j]=0.5*(pow(std::abs(M(0,0)),2.0)+pow(std::abs(M(0,1)),2.0)
                +pow(std::abs(M(1,0)),2.0)+pow(std::abs(M(1,1)),2.0));
        }
}

void ExtendedJones::calculateOneLambdaWithStokes(int iLambda){
    //If the program comes here, it means there is at least one polarizer and its theta angle approaches to 90 or 270 degree.

}

const TRANSRESULT& ExtendedJones::getTransmissions(){
    return transmissions;
}

const STOKESRESULT& ExtendedJones::getStokes(){
    return stokes;
}

void ExtendedJones::resetTransmissions(){
    transmissions.clear();
    transmissions = DOUBLEARRAY2D(_inAngles.size(), DOUBLEARRAY1D(_inAngles[0].size(), 0.0));
}
void ExtendedJones::resetTransTemp(){
    transTemp.clear();
    transTemp = DOUBLEARRAY2D(_inAngles.size(), DOUBLEARRAY1D(_inAngles[0].size(), 0.0));
}
void ExtendedJones::resetStokes(){
    stokes.clear();
    stokes = STOKESRESULT(lambdas.size(), std::vector<std::vector<STOKESTRACE> >(_inAngles.size(),
        std::vector<STOKESTRACE>(_inAngles[0].size())));
}

void ExtendedJones::resetToCalculateWithNewDiretors(DIRVEC _in){
    if (lcLayerindex >= 0){
        std::shared_ptr<Optical2X2OneLayer<UniaxialType> > tempLayerPtr;
        tempLayerPtr = std::dynamic_pointer_cast<Optical2X2OneLayer<UniaxialType>,
            Optical2X2OneLayerBase>(matLayers[lcLayerindex]);
        if (tempLayerPtr) tempLayerPtr->resetDirectors(_in);
    }
}
