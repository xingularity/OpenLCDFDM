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

#ifndef LCD1D_LAYERFACTORY_HPP
#define LCD1D_LAYERFACTORY_HPP

#include "LCD_Optics2x2.hpp"

namespace LCD1D{
    using LCDOptics::OpticalMaterialClass;
    using LCD::COMPD;
    using LCDOptics::Optical2x2BasePtr;
    using LCDOptics::Optical2x2IsoPtr;
    using LCDOptics::Optical2x2UnixialPtr;
     ///This class is used to generate all optical layers.
    class OpticalLayerFactory{
    public:
        OpticalLayerFactory();
        void generateIsotropicLayer(double thickness, std::map<double, std::vector<COMPD> > dispersion, OpticalMaterialClass layerKind = OpticalMaterialClasss::ISOTROPIC);
        void generateUniaxialLayer(double thickness, std::map<double, std::vector<COMPD> > dispersion, DIRVEC axes,
            OpticalMaterialClass layerKind = OpticalMaterialClasss::UNIAXIAL);
        void generateBiaxialLayer(double thickness, std::map<double, std::vector<COMPD> > dispersion, DIRVEC axes,
            OpticalMaterialClass layerKind = OpticalMaterialClasss::BIAXIAL);
        MATERIALLAYERS2X2CONT getOpticalLayers();
    private:
        MATERIALLAYERS2X2CONT matLayers;
    };
};

#endif
