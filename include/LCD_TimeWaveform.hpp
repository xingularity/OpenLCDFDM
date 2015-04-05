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
 *   * Neither the name of OpenLCDFDM nor
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

#ifndef LCD_TIMEWAVEFORM_HPP
#define LCD_TIMEWAVEFORM_HPP

#include "LCD_ContainerDefine.hpp"
#include <map>

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
        StepWaveform(std::map<double, double> _prof, double _period, double _ampShift=0.0, double _preApplyTime = 0.0, double _preApplyVal=0.0);
        virtual double operator()(double t);
        DOUBLEARRAY2D getStepProfile();
        protected:
        double wavePeriod{0.0};
        int stepNum{0};
        DOUBLEARRAY2D prof;
        double ampShift{0.0};
    };
};

#endif
