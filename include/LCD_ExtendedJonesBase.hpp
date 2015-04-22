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

#ifndef LCD_EXTENDEDJONESBASE_HPP
#define LCD_EXTENDEDJONESBASE_HPP

#include "LCD_Optics2x2.hpp"
#include <cassert>
#include <limits>

namespace LCDOptics{
    using LCD::DOUBLEARRAY2D;
    using Optical2x2BasePtr = std::shared_ptr<Optical2X2OneLayerBase>;
    using Optical2x2IsoPtr = std::shared_ptr<Optical2X2OneLayer<ISOType> >;
    using Optical2x2UnixialPtr = std::shared_ptr<Optical2X2OneLayer<UniaxialType> >;
    ///container of all optical materials, in the same order of passing through a ray
    typedef std::vector<std::shared_ptr<Optical2X2OneLayerBase> > MATERIALLAYERS2X2CONT;
    ///spectral efficiency of human eyes, get it from yBar of lambda from CIE web site.
    const LIGHTSPECTRUMDATA SpectrumEfficiency={
        {0.38, 0.000039},
        {0.385, 0.000064},
        {0.39, 0.00012},
        {0.395, 0.000217},
        {0.400, 0.000396},
        {0.405, 0.00064},
        {0.410, 0.00121},
        {0.415, 0.00218},
        {0.420, 0.004},
        {0.425, 0.0073},
        {0.430, 0.0116},
        {0.435, 0.01684},
        {0.440, 0.023},
        {0.445, 0.0298},
        {0.450, 0.038},
        {0.455, 0.048},
        {0.460, 0.06},
        {0.465, 0.0739},
        {0.470, 0.09098},
        {0.475, 0.1126},
        {0.480, 0.13902},
        {0.485, 0.1693},
        {0.490, 0.20802},
        {0.495, 0.2586},
        {0.500, 0.323},
        {0.505, 0.4073},
        {0.510, 0.503},
        {0.515, 0.6082},
        {0.520, 0.71},
        {0.525, 0.7932},
        {0.530, 0.862},
        {0.535, 0.91485},
        {0.540, 0.954},
        {0.545, 0.9803},
        {0.550, 0.99495},
        {0.555, 1},
        {0.560, 0.995},
        {0.565, 0.9786},
        {0.570, 0.952},
        {0.575, 0.9154},
        {0.580, 0.87},
        {0.585, 0.8163},
        {0.590, 0.757},
        {0.595, 0.6949},
        {0.600, 0.631},
        {0.605, 0.5668},
        {0.610, 0.503},
        {0.615, 0.4412},
        {0.620, 0.381},
        {0.625, 0.321},
        {0.630, 0.265},
        {0.635, 0.217},
        {0.640, 0.175},
        {0.645, 0.1382},
        {0.650, 0.107},
        {0.655, 0.0816},
        {0.660, 0.061},
        {0.665, 0.04458},
        {0.670, 0.032},
        {0.675, 0.0232},
        {0.680, 0.017},
        {0.685, 0.01192},
        {0.690, 0.00821},
        {0.695, 0.005723},
        {0.700, 0.004102},
        {0.705, 0.002929},
        {0.710, 0.002091},
        {0.715, 0.001484},
        {0.720, 0.001047},
        {0.725, 0.00074},
        {0.730, 0.00052},
        {0.735, 0.000361},
        {0.740, 0.000249},
        {0.745, 0.000172},
        {0.750, 0.00012},
        {0.755, 0.000085},
        {0.760, 0.00006},
        {0.765, 0.000042},
        {0.770, 0.00003},
        {0.775, 0.000021},
        {0.780, 0.000015}
    };

    /**
    Base class of the main loop to calculate extended Jones matrix.
    */
    class ExtendedJonesBase{
    public:
        ///For sigle wavelength calculation, one need no LIGHTSPECTRUMDATA in single wavelength, if input, it will be ignored.
        ExtendedJonesBase(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles,const double targetLambda, 
            LIGHTSPECTRUMDATA lightSrcSpectrum_ = LIGHTSPECTRUMDATA()):
        matLayers(_materials), inAngles(_inAngles), lambdas(1, targetLambda)
        {
            matLayerNum = matLayers.size();
            //Don't actually need these for single wave length calculation
            lightSourceSpectrum = DOUBLEARRAY1D(1, 1.0);
            yBarOfLambda = DOUBLEARRAY1D(1, 1.0);
            if (matLayerNum < 1){
                std::cout << "Can't calculate optics without an layer of optical material." << std::endl;
                assert(false);
            }
            checkPolazerLayers();
        }
        ///For multiwavelengths calculation
        ExtendedJonesBase(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
        const double end_lambda_, const double step_lambda_, LIGHTSPECTRUMDATA lightSrcSpectrum_ = LIGHTSPECTRUMDATA())
        :matLayers(_materials), inAngles(_inAngles){
            //create lambdas first
            lambdas.clear();
            for (double ll = start_lambda_; ll <= end_lambda_ + 10.0*std::numeric_limits<double>::epsilon(); ll += step_lambda_)
                lambdas.push_back(ll);
            if (lambdas.size() < 2){
                std::cout << "lambdas.size() < 2 for multi-wavelength calculation." << std::endl;
                assert(false);
            }
            matLayerNum = matLayers.size();
            if (matLayerNum < 1){
                std::cout << "Can't calculate optics without any layer of optical material." << std::endl;
                assert(false);
            }
            // if input lightSrcSpectrum_ is empty, then use equal energy spectrum
            if (lightSrcSpectrum_.size() == 0){
                lightSrcSpectrum_[0.38] = 1.0;
                lightSrcSpectrum_[0.78] = 1.0;
            }
            //begin to do spectrum interpolation
            SpectrumInterpolator<LIGHTSPECTRUMDATA> interpolator(lambdas);
            interpolator.interpolate(SpectrumEfficiency, yBarOfLambda);
            interpolator.interpolate(lightSrcSpectrum_, lightSourceSpectrum);
            checkPolazerLayers();
        }

        const int getLCLayerIndex(){return lcLayerindex;}
        const MATERIALLAYERS2X2CONT& getMaterialLayers()const {return matLayers;}
        const IAngles& getIncidentAngles()const{return inAngles;}
        const DOUBLEARRAY1D& targetLambdas()const{return lambdas;}
        const DOUBLEARRAY1D& targetLightSrcSpectrum()const{return lightSourceSpectrum;}
        const DOUBLEARRAY1D& targetYBarOfLambda()const{return yBarOfLambda;}

    protected:
        void findLCLayerInMaterialList(){
            //find LC layer
            for (int i = 0; i < matLayers.size(); ++i){
                Optical2X2OneLayer<UniaxialType>* sp = dynamic_cast<Optical2X2OneLayer<UniaxialType>*> (matLayers[i].get());
                if (sp == NULL) continue;
                if (sp->opticalLayerKind() == OPT_LCMATERIAL) {
                    lcLayerindex = i;
                    break;
                }
            }
        }
         ///Stokes calculation needs to know the positions of polarizer layer.
        void checkPolazerLayers(){
            for (int i = 0; i < matLayers.size(); ++i)
                if (matLayers[i]->opticalLayerKind() == OPT_POLARIZER){
                    polarizerLayersIndex.push_back(i);
                }
        }
        ///Incident angles which will be calculated.
        const IAngles inAngles;
        ///material list
        MATERIALLAYERS2X2CONT matLayers;
        ///size of materials
        int matLayerNum;
        ///calculation wavelengths
        DOUBLEARRAY1D lambdas;
        ///Spectrum of light source, interpolated to target wavelengths
        DOUBLEARRAY1D lightSourceSpectrum;
        ///spectral efficiency of human eyes, interpolated to target wavelengths.
        DOUBLEARRAY1D yBarOfLambda;
        ///refractive index of the air
        const double nAir{1.000293};
        int lcLayerindex{-1};
        std::vector<unsigned int> polarizerLayersIndex;
    };
};

#endif
