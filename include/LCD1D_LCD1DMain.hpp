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

#ifndef LCD1D_LCD1DMAIN_HPP
#define LCD1D_LCD1DMAIN_HPP

#include <omp.h>
#include <ctime>
#include <stdexcept>
#include "LCD1D_ExtendedJones.hpp"
#include "LCD_UsefulFuncs.hpp"
#include "LCD1D_FDM1DSolver.hpp"

class LCD1DMainBase{
public:
    LCD1DMainBase(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing);
    void setTFTPI(LCD1D::DielecParameters _tftpi);
    void setCFPI(LCD1D::DielecParameters _cfpi);
    void setOMPThreadNum(size_t _num);
    
protected:
    std::share_ptr<LCD1D::Epsilon> epsilonr;
    std::share_ptr<LCD1D::Potential> potentials;
    std::share_ptr<LCD1D::LCDirector> lcDir;
};

class LCD1DStaticMain: public LCD1DMainBase{
public:
    LCD1DStaticMain(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing, 
        double _voltStart, double _voltEnd, double _voltStep);
    
};

class LCD1DDynamicMain: public LCD1DMainBase{
public:
    LCD1DStaticMain(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing, double _maxCalcTime);
    void setDCWaveform(double _volt);
    void setStepWaveform(std::map<double, double> _profile);
    
private:
    std::shared_ptr<WaveformBase> waveform;
    double maxCalcTime;
};

#endif
