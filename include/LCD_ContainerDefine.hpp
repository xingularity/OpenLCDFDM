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
 *   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
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
