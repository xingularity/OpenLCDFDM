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

#ifndef LCD1D_FDM1DSOLVER_HPP
#define LCD1D_FDM1DSOLVER_HPP

#include "LCD_ContainerDefine.hpp"

namespace LCD1D{
    ///Rubbing conditions are tft-side theta, tft-side phi, cf-side theta, total twist.
    struct LCConditions{
        double tftTheta;
        double tftPhi;
        double cfTheta;
        double totalTwist;
        double epsPara;
        double epsPerp;
        double gamma;
        double k11;
        double k22;
        double k33;
    };
    
    class LCDirector{
        public:
            ///Constructor 
            LCDirector(const LCConditions cond_, size_t size_);
            const getSize()const{return size_;}
            const DIRVEC getDirectors()const {return dir;}
        private:
            LCD::DIRVEC dir;
    };
    
    class Potential{
        public:
            Potential(int size);
        private:
            
    };
};

#endif
