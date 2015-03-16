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

#include "LCD_TimeWaveform.hpp"

DCWaveform::DCWaveform(double val, double _preApplyTime, double _preApplyVal){
    dc = val;
    preApplyTime = _preApplyTime;
    preApplyVal = _preApplyVal;
}

double DCWaveform::operator()(double t){
    return dc;
}

StepWaveform::StepWaveform(std::map<double, double> _prof, double _period, double _ampShift, 
                           double _preApplyTime, double _preApplyVal){
    period = _period;
    stepNum = _prof.size();
    wavePeriod = _period;
    ampShift = _ampShift;
    preApplyTime = _preApplyTime;
    preApplyVal = _preApplyVal;
    prof = DOUBLEARRAY2D(stepNum, DOUBLEARR1D(2, 0.0));
    int index = 0;
    for (auto i : _prof){
        prof[index][0] = i.first;
        prof[index][1] = i.second;
        ++index;
    }
    //if the given profile doesn't start with 0 time, then shift it.
    if (std::abs(prof[0][0] - 0.0) >= 1.0e-13){
        int shift = prof[0][0];
        for (auto i : prof){
            i[0] -= shift;
        }
    }
}

double StepWaveform::operator()(double t){
    if (t < preApplyTime){return preApplyVal;}
    
    //t <= preApplyTime isn't regarded as part of the step waveform
    t -= preApplyTime;
    
    int temp=0;

    temp=std::floor((t-prof[0][0])/period);
    if ((wavePeriod>1.0e-13) && (temp >= 1.0)){t=t-temp*stepPeriod;}

    if (t >= prof[stepNum-1][0]){
        return prof[stepNum-1][1] + ampShift;
    }

    for (int i=0; i< stepNum-1; i++){
        if ((t-prof[i][0])*(t-prof[i+1][0]) <= 1.0e-13){
            return prof[i][1]+ampShift;
        }
    }
    return 0.0;

}
