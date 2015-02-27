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
    using LCDOptics::MATERIALLAYERS2X2CONT;
    using LCDOptics::IAngles;
    using LCDOptics::Angle;
    using LCDOptics::LIGHTSPECTRUMDATA;
    using LCD::DIRVEC;
    using LCDOptics::JONESMAT;
    using LCD::DOUBLEARRAY1D;
    using LCD::DOUBLEARRAY2D;
    using TRANSRESULT = DOUBLEARRAY2D;
    using STOKESTRACE = std::vector<Eigen::Vector3d>;
    using STOKESRESULT = std::vector<std::vector<std::vector<STOKESTRACE> > >;
    ///Main function for extended Jones matrix
    class ExtendedJones: public LCDOptics::ExtendedJonesBase{
    public:
        ///For sigle wavelength calculation
        ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles,const double targetLambda,
            LIGHTSPECTRUMDATA lightSrcSpectrum_, bool _ifLambertian = false, bool _ifStokes=false);
        ///For multiwavelengths calculation
        ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
        const double end_lambda_, const double step_lambda_, LIGHTSPECTRUMDATA lightSrcSpectrum_,
        bool _ifLambertian = false, bool _ifStokes=false);
        ///main function to calculate extended Jones.
        void calculateExtendedJones();
        ///return the transmission results
        const DOUBLEARRAY2D& getTransmissions();
        ///return the results of stokes values
        const STOKESRESULT& getStokes();
        ///reset directors in materials and the states of computing components to calculate optics based on input directors.
        void resetLCDiretors(DIRVEC _in);
        ///return if this object calculates stokes values
        bool isCalcStokes()const {return ifCalcStokes;}
        ///set number of threads for OpenMP
        void setNumThreads(int _numThreads);
    private:
        ///If LC layer is sandwiched between 2 polarizer layers with the theta angles of their optical axis equal to 90 or 270 degree, the stokes value can be calculated.
        void checkIfCalcStokes();
        ///for single-wavelength calculation, no stokes calculation.
        void calculateOneLambdaNoStokes(int iLambda);
        ///for single-wavelength calculation and stokes calculation.
        void calculateOneLambdaWithStokes(int iLambda);
        ///for multi-wavelength calculation
        void calculateManyLambda();
        void resetTransmissions();
        void resetTransEachLambda();
        void resetStokes();
        ///results of transmissions on corresponding angles in inAngles
        DOUBLEARRAY2D transmissions;
        ///array of transmission results of every wavelength
        std::vector<DOUBLEARRAY2D> transEachLambda;
        ///store the stokes values
        STOKESRESULT stokes;
        ///if user wanna calculate Stokes
        bool ifCalcStokes{false};
        bool ifLambertian{false};
        ///
        unsigned int ompNumThreads{0};
    };
};

#endif
