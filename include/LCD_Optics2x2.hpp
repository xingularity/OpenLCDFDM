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
        ///nk dispersion for isotropic material
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
        void resetDirectors(DIRVEC _in);
        ///calculate one jones matrix and return its average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, Angle& iang, double lambda, double lastn);
        ///calculate one jones matrix and light polarization, return the average refractive index and incident angles
        virtual double calcJonesMatrix(JONESMAT& m, POLARTRACE& lightPolar, Angle& iang, double lambda, double lastn);
        virtual void interpolateNKForLambdas(DOUBLEARRAY1D lambdas_);
    private:
        ///find nk for corresponding lambda. If find in map, then return. Otherwise, interpolate immediately
        NKoNKe findNK(double lambda);
        ///nk dispersion for uniaxial material
        NKoNKeData nk;
        ///Interpolation of nk distribution
        SpectrumInterpolator<NKoNKeData> nkInterpolator;
    };
};

#endif
