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

#ifndef LCD1D_EXTENDEDJONES_HPP
#define LCD1D_EXTENDEDJONES_HPP

#include "LCD_ExtendedJonesBase.hpp"

namespace LCD1D{
    using LCDOptics::POLARTRACE;
    using TRANSRESULT = DOUBLEARRAY2D;
    using STOKESTRACE = std::vector<Eigen::Vector3d>;
    using STOKESRESULT = std::vector<std::vector<std::vector<STOKESTRACE> > >;
    ///Main function for extended Jones matrix
    class ExtendedJones: LCDOptics::ExtendedJonesBase{
    public:
        ///For sigle wavelength calculation
        ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles,const double targetLambda,
            LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifStokes=false);
        ///For multiwavelengths calculation
        ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
        const double end_lambda_, const double step_lambda_, LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifStokes=false);
        ///main function to calculate extended Jones.
        void calculateExtendedJones();
        ///return the transmission results
        const TRANSRESULT& getTransmissions();
        ///return the results of stokes values
        const STOKESRESULT& getStokes();
        ///reset directors in materials and the states of computing components to calculate optics based on input directors.
        void resetToCalculateWithNewDiretors(DIRVEC _in);
        ///return if this object calculates stokes values
        bool isCalcStokes()const {return ifCalcStokes;}
    private:
        ///for single-wavelength calculation, no stokes calculation.
        void calculateOneLambdaNoStokes(int iLambda);
        ///for single-wavelength calculation and stokes calculation.
        void calculateOneLambdaWithStokes(int iLambda);
        ///for multi-wavelength calculation
        void calculateManyLambda();
        void resetTransmissions();
        void resetTransTemp();
        void resetStokes();
        ///results of transmissions on corresponding angles in inAngles
        TRANSRESULT transmissions;
        ///temporary array for transmission results of one wave length for multi-wavelength calculation.
        TRANSRESULT transTemp;
        ///store the stokes values
        STOKESRESULT stokes;
        ///if user wanna calculate Stokes
        bool ifCalcStokes{false};
    };
};

#endif
