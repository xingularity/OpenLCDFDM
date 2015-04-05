/*
 * Copyright (C) 2015 Zong-han, Xie <icbm0926@gmail.com>.
 *
 * You may use this file under the terms of the BSD license as follows:
 *
 * "Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of OpenLCDFDM nor
 *     the names of its contributors may be used to endorse or promote
 *     products derived from this software without specific prior written
 *     permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
 *
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
    };
};

#endif
