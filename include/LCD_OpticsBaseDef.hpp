/*
 * Copyright (C) 2015 Zong-han, Xie <icbm0926@gmail.com>.
 *
 * You may use this file under the terms of the BSD license as follows:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of OpenLCDFDM nor the names of its contributors 
 *     may be used to endorse or promote products derived from this 
 *     software without specific prior written permission.
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
 
#ifndef LCD_OPTICSBASEDEF_HPP
#define LCD_OPTICSBASEDEF_HPP

#include "LCD_ContainerDefine.hpp"
#include <memory>
#include <map>
#include <tuple>
#include <stdexcept>

namespace LCDOptics{
    using LCD::BZARRAY;
    using LCD::DOUBLEARRAY1D;
    using LCD::COMPLEXDARRAY1D;
    using LCD::COMPD;
    using LCD::DIRVEC;
    using std::runtime_error;
    ///incident angle data, (theta, phi)
    using Angle = std::tuple<double, double>;
    ///incident angles
    using IAngles = std::vector<std::vector<Angle> >;
    ///record a set of No, Ne data
    using NKoNKe = std::tuple<COMPD,COMPD>;
    ///record a set of N1, N2, N3 data
    using NK123 = std::tuple<COMPD,COMPD,COMPD>;
    ///type of the input spectrum data of incident light
    using LIGHTSPECTRUMDATA = std::map<double, double>;
    ///type of input nk dispersion
    using NKData = std::map<double, COMPD>;
    ///type of input nk dispersion
    using NKoNKeData = std::map<double,NKoNKe>;
    ///type of input nk dispersion
    using NK123Data = std::map<double,NK123>;
    ///enum for optical material type
    enum OpticMaterialType{
        ///Isotropic optical material layer
        ISOType=0,
        ///Uniaxia optical material layer
        UniaxialType=1,
        ///Biaxial optical material layer
        BiaxialType=2,
    };

    ///create an "Angle" tuple with radian data
    inline Angle makeAngle(double theta, double phi){
        return std::make_tuple(theta, phi);
    }

    ///create an "Angle" tuple with degree data
    inline Angle makeAngle2(double theta, double phi){
        return std::make_tuple(theta*M_PI/180.0, phi*M_PI/180.0);
    }

    /**
    SpectrumInterpolator interpolates given nk dispersion to taget nk dispersion.<br/>
    Using linear interpolation.
    NKINPUTDATA must be a std::map, it can be one of LCDOptics::NKData, LCDOptics::NKoNKeData, LCDOptics::NK123Data, LCDOptics::LIGHTSPECTRUMDATA.
    */
    template <typename INTYPE>
    class SpectrumInterpolator{
    public:
        SpectrumInterpolator();
        /// constructor, input target wavelengths.
        SpectrumInterpolator(const DOUBLEARRAY1D lambda_);
        ///interpolate input nk to the distribution on target
        void interpolate(const INTYPE& nkData, std::vector<typename INTYPE::value_type::second_type>& wanted_nk);
        void interpolate(const INTYPE& nkData, INTYPE& wanted_nk);
        ///only interpolate for one lambda
        typename INTYPE::value_type::second_type interpolate(const INTYPE& nkData, double lambda);
        ///reset target wavelengths
        void resetTargetLambdas(const DOUBLEARRAY1D lambda_){lambdas = lambda_;};
        DOUBLEARRAY1D& getLambdas(){return lambdas;}
    private:
        ///target wavelengths to interpolate
        DOUBLEARRAY1D lambdas;
        /**
        Find the range where the wavelength lies within.
        If given lambda is smaller than the smallest value in map, the returned tuple contains two same values, they all point to the smallest key in map.
        If given lambda is larger than the largest value in map, the returned tuple contains two same values, they all point to the largest key in map.
        */
        std::tuple<double, double> searchLambdaRegion(const double lambda, const INTYPE& in);
        ///toleratedError is a value specifying tolerated floating point error for comparing wavelegnth values while looking for the range in map.
        double toleratedError;
    };
};
#endif
