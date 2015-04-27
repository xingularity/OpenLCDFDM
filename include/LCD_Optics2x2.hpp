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
        void resetAxes(DIRVEC _in){
            layernum = _in.size();
            axisVec.resize(_in.size());
            axisVec = _in;
        }
        void resetThickness(double _thick){d = thick;}
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
    class Optical2X2OneLayer: public Optical2X2OneLayerBase{
        Optical2X2OneLayer(){throw runtime_error("Can't use this Optical2X2OneLayer");}
    };

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
        ///input nk dispersion for isotropic material
        const NKData nk_in;
        ///nk dispersion used in calculation
        NKData nk;
        ///Interpolation of nk distribution
        SpectrumInterpolator<NKData> nkInterpolator;
    };

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
        ///calculate one jones matrix and return its average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, Angle& iang, double lambda, double lastn);
        ///calculate one jones matrix and light polarization, return the average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, POLARTRACE& lightPolar, Angle& iang, double lambda, double lastn);
        virtual void interpolateNKForLambdas(DOUBLEARRAY1D lambdas_);
    private:
        ///find nk for corresponding lambda. If find in map, then return. Otherwise, interpolate immediately
        NKoNKe findNK(double lambda);
        ///input nk dispersion for uniaxial material
        const NKoNKeData nk_in;
        ///nk dispersion for calculation
        NKoNKeData nk;
        ///Interpolation of nk distribution
        SpectrumInterpolator<NKoNKeData> nkInterpolator;
    };
};

#endif
