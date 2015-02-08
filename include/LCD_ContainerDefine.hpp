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

#ifndef LCD_CONTAINERDEFINE_HPP
#define LCD_CONTAINERDEFINE_HPP

#include "blitz/array.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <complex>
#include <cmath>

namespace LCD{
    template<typename T, int N>
    using BZARRAY = blitz::Array<T, N>;
    template<typename T, int N>
    using BZTINYVEC = blitz::TinyVector<T, N>;

    using COMPD = std::complex<double>;
    using EigenM22 = Eigen::Matrix2d;
    using EigenM44 = Eigen::Matrix4d;
    using EigenC22 = Eigen::Matrix2cd;
    using EigenC44 = Eigen::Matrix4cd;
    ///type to contain a director or a 3D vector using double
    using BzVECD3D = BZTINYVEC<double, 3> ;
    ///type to contain an array of director or a 3D vector using double
    using DIRVEC = BZARRAY<BzVECD3D, 1> ;
    ///type to contain an array of director or a 3D vector using double
    using DIRMATRIX2D = BZARRAY<BzVECD3D, 2> ;
    ///type to contain an array of director or a 3D vector using double
    using DIRMATRIX3D = BZARRAY<BzVECD3D, 3> ;

    using INTARRAY1D = std::vector<int> ;
    using DOUBLEARRAY1D = std::vector<double> ;
    using DOUBLEARRAY2D = std::vector<std::vector<double> >;
    using LONGARRAY1D = std::vector<long> ;
    using COMPLEXDARRAY1D = std::vector<COMPD> ;

};

#endif
