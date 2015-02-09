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

#ifndef LCD1D_EXTENDEDJONES1D_HPP
#define LCD1D_EXTENDEDJONES1D_HPP

#include "LCD_ExtendedJonesBase.hpp"

namespace LCD1D{
    using TRANSRESULT = DOUBLEARRAY2D;
    using STOKESRESULT = std::vector<std::vector<std::vector<Eigen::Vector3d> > >;
    class ExtendedJones: LCDOptics::ExtendedJonesBase{
    public:
        ///For sigle wavelength calculation
        ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles,const double targetLambda, LIGHTSPECTRUMDATA lightSrcSpectrum_);
        ///For multiwavelengths calculation
        ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
        const double end_lambda_, const double step_lambda_, LIGHTSPECTRUMDATA lightSrcSpectrum_);
        void calculateExtendedJones(bool ifStokes = false);
        const TRANSRESULT& getTransmissions();
        const STOKESRESULT& getStokes();
    private:
        ///results of transmissions on corresponding angles in inAngles
        TRANSRESULT transmissions;
        ///temporary array for transmission results of one wave length for multi-wavelength calculation.
        TRANSRESULT transTemp;
        ///store the stokes values
        STOKESRESULT stokes;
    };
};

#endif
