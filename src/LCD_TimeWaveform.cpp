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

#include "LCD_TimeWaveform.hpp"

using namespace LCD;

DCWaveform::DCWaveform(double val, double _preApplyTime, double _preApplyVal){
    dc = val;
    preApplyTime = _preApplyTime;
    preApplyVal = _preApplyVal;
}

double DCWaveform::operator()(double t){
    if (t <= preApplyTime)
        return preApplyVal;
    return dc;
}

StepWaveform::StepWaveform(std::map<double, double> _prof, double _period, double _ampShift,
                           double _preApplyTime, double _preApplyVal){
    wavePeriod = _period;
    stepNum = _prof.size();
    wavePeriod = _period;
    ampShift = _ampShift;
    preApplyTime = _preApplyTime;
    preApplyVal = _preApplyVal;
    prof = DOUBLEARRAY2D(stepNum, DOUBLEARRAY1D(2, 0.0));
    size_t index = 0;
    for (auto i : _prof){
        prof[index][0] = i.first;
        prof[index][1] = i.second;
        ++index;
    }
    //if the given profile doesn't start with 0 time, then shift it.
    if (std::abs(prof[0][0] - 0.0) >= 1.0e-13){
        double shift = prof[0][0];
        for (auto& i : prof)
            i[0] -= shift;
    }
}

double StepWaveform::operator()(double t){
    if(wavePeriod < 1.0e-13) {return 0.0;}

    if (t < preApplyTime){return preApplyVal;}

    //t <= preApplyTime isn't regarded as part of the step waveform
    t -= preApplyTime;

    int temp=0;

    temp=std::floor((t-prof[0][0])/wavePeriod);
    if ((wavePeriod>1.0e-13) && (temp >= 1.0)){t=t-temp*wavePeriod;}

    if (t >= prof.back()[0]){
        return prof.back()[1] + ampShift;
    }

    for (int i=0; i< stepNum-1; i++){
        if (std::abs(t-prof[i][0]) <= 1.0e-13) return prof[i][1]+ampShift;
        if (std::abs(t-prof[i+1][0]) <= 1.0e-13) return prof[i+1][1]+ampShift;
        if ((t-prof[i][0])*(t-prof[i+1][0]) <= 1.0e-13)return prof[i][1]+ampShift;
    }
}

DOUBLEARRAY2D StepWaveform::getStepProfile(){
    return prof;
}
