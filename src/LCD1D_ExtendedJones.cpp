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
#include "LCD1D_ExtendedJones.hpp"

using namespace LCD1D;

ExtendedJones::ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles,const double targetLambda, std::map<double, double> lightSrcSpectrum_):
ExtendedJonesBase(_materials, _inAngles, targetLambda, lightSrcSpectrum_)
{
}

ExtendedJones::ExtendedJones(MATERIALLAYERS2X2CONT& _materials, const IAngles _inAngles, const double start_lambda_,
const double end_lambda_, const double step_lambda_, std::map<double, double> lightSrcSpectrum_)
ExtendedJonesBase(_materials, _inAngles, start_lambda_, end_lambda_, step_lambda_, lightSrcSpectrum_)
{
}