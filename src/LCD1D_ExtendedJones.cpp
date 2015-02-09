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
    resetTransmissions();
    resetTransTemp();
    findLCLayerInMaterilList();
    ifCalcStokes = _ifStokes;
    if (ifCalcStokes)resetStokes();
}

ExtendedJones::ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
const double end_lambda_, const double step_lambda_, LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifStokes)
ExtendedJonesBase(_materials, _inAngles, start_lambda_, end_lambda_, step_lambda_, lightSrcSpectrum_)
{
    resetTransmissions();
    resetTransTemp();
    findLCLayerInMaterilList();
    ifCalcStokes = _ifStokes;
    if (ifCalcStokes)resetStokes();
}

void ExtendedJones::calculateExtendedJones(bool ifStokes = false){

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
    if (lambdas.size() < 2) return;
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
    resetTransmissions();
    resetTransTemp();
    resetStokes();
}
