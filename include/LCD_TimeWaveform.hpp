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

#ifndef LCD_TIMEWAVEFORM_HPP
#define LCD_TIMEWAVEFORM_HPP

#include "blitz/array.h"
#include "LCD_ContainerDefine.hpp"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <complex>
#include <cmath>

namespace LCD{
    enum WaveformType{
        Waveform_DCType=0,
        Waveform_StepType=1,
    };

    class WaveformBase{
        public:
        WaveformBase(){};
        virtual double operator()(double t) = 0;
        WaveformType type(){return waveType;}
        protected:
        double preApplyTime{0.0};
        double preApplyVal{0.0};
        WaveformType waveType;
    };
    
    class DCWaveform:public WaveformBase{
        public:
        DCWaveform(double val, double _preApplyTime=0.0, double _preApplyVal=0.0);
        virtual double operator()(double t);
        protected:
        double dc{0.0};
    };

    class StepWaveform:public WaveformBase{
        public:
        StepWaveform(std::map<double, double> _prof, double _period, double _ampShift, double _preApplyTime = 0.0, double _preApplyVal=0.0);
        virtual double operator()(double t);
        protected:
        double wavePeriod{0.0};
        int stepNum{0};
        DOUBLEARRAY2D prof;
        double ampShift{0.0};
    };
    
};

#endif
