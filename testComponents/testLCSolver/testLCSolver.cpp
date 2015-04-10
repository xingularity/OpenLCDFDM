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

#include "LCD1D_FDM1DSolver.hpp"
#include <limits>
#include <fstream>
#include <cmath>
#include <set>
#include "LCD_UsefulFuncs.hpp"

void writeDirectors(LCD1D::DIRVEC dirs, std::string suffix=""){
    std::fstream file;
    file.open(("Directors"+suffix+".txt").c_str(), std::ios::out|std::ios::trunc);
    for (int i = 0; i < dirs.extent(0); ++i){
        file << dirs(i)(0) << ", " << dirs(i)(1) << ", " << dirs(i)(2) << ", " << std::acos(dirs(i)(2))
             << ", " << std::atan2(dirs(i)(1), dirs(i)(0)) << " " << std::endl;
    }
    file.flush();
    file.close();
}

void writePotentials(LCD::DOUBLEARRAY1D potentials, std::string suffix=""){
    std::fstream file;
    file.open(("Potentials"+suffix+".txt").c_str(), std::ios::out|std::ios::trunc);
    for (int i = 0; i < potentials.size(); ++i){
        file << potentials[i] << std::endl;
    }
    file.flush();
    file.close();
}

void writeEFieldForLC(LCD::DOUBLEARRAY1D EFieldForLC, std::string suffix=""){
    std::fstream file;
    file.open(("EFieldForLC"+suffix+".txt").c_str(), std::ios::out|std::ios::trunc);
    for (int i = 0; i < EFieldForLC.size(); ++i){
        file << EFieldForLC[i] << std::endl;
    }
    file.flush();
    file.close();
}

void testTNCalculation(double dc_volt){
    std::cout << "Start TN Calculation without q0 under dc_volt=" + toString(dc_volt) << std::endl;
    size_t layerNum = 40;
    double dt = 0.01; //ms
    size_t maxIter = 100000000;
    double max_error = 1.0e-8;
    LCD1D::LCParamters lcParam = {4.0, 12, 3.6, 60, 12.0, 6.5, 15.0, 0.0};
    LCD1D::DielecParameters tftpiParam = {0.1, 3.5};
    LCD1D::DielecParameters cfpiParam = {0.1, 3.5};
    LCD1D::RubbingCondition rubbing = {89.0*M_PI/180.0, 45.0*M_PI/180.0, 89.0*M_PI/180.0, 90.0*M_PI/180.0};
    LCD1D::Epsilon epsilonr(layerNum, lcParam.thick);
    epsilonr.setTFTPI(tftpiParam);
    epsilonr.setCFPI(cfpiParam);
    LCD1D::Potential potential(layerNum+1, epsilonr, lcParam.thick/layerNum);
    potential.createStaticUpdatePolicy();
    LCD1D::LCDirector lcDirector(layerNum, lcParam, rubbing, epsilonr);
    lcDirector.createVectorFormUpdater(potential, dt);
    double residual = std::numeric_limits<double>::max();
    //writeDirectors(lcDirector.getDirectors(), "_TNNoq0_"+toString(dc_volt)+"V_initialized_");
    while(residual >= max_error){
        potential.update(dc_volt);
        residual = lcDirector.update();
        std::cout << "residual = " << residual << std::endl;
    };
    writePotentials(potential.getPotentials(), "_TNNoq0_"+toString(dc_volt)+"V");
    //writeEFieldForLC(potential.getEfieldForLC(), "_TNNoq0_"+toString(dc_volt)+"V");
    writeDirectors(lcDirector.getDirectors(), "_TNNoq0_"+toString(dc_volt)+"V");
}

void testTNWithQ0Calculation(double dc_volt){
    std::cout << "Start TN Calculation with q0 under dc_volt=" + toString(dc_volt) << std::endl;
    size_t layerNum = 40;
    double dt = 0.01; //ms
    size_t maxIter = 100000000;
    double max_error = 1.0e-8;
    LCD1D::LCParamters lcParam = {4.0, 12, 3.6, 60, 12.0, 6.5, 15.0,  2.0*M_PI/70.0};
    LCD1D::DielecParameters tftpiParam = {0.1, 3.5};
    LCD1D::DielecParameters cfpiParam = {0.1, 3.5};
    LCD1D::RubbingCondition rubbing = {89.0*M_PI/180.0, 45.0*M_PI/180.0, 89.0*M_PI/180.0, 90.0*M_PI/180.0};
    LCD1D::Epsilon epsilonr(layerNum, lcParam.thick);
    epsilonr.setTFTPI(tftpiParam);
    epsilonr.setCFPI(cfpiParam);
    LCD1D::Potential potential(layerNum+1, epsilonr, lcParam.thick/layerNum);
    potential.createStaticUpdatePolicy();
    LCD1D::LCDirector lcDirector(layerNum, lcParam, rubbing, epsilonr);
    lcDirector.createVectorFormUpdater(potential, dt);
    double residual = std::numeric_limits<double>::max();
    //writeDirectors(lcDirector.getDirectors(), "_TNNoq0_"+toString(dc_volt)+"V_initialized_");
    while(residual >= max_error){
        potential.update(dc_volt);
        residual = lcDirector.update();
        std::cout << "residual = " << residual << std::endl;
    };
    writePotentials(potential.getPotentials(), "_TNq0_"+toString(dc_volt)+"V");
    //writeEFieldForLC(potential.getEfieldForLC(), "_TNNoq0_"+toString(dc_volt)+"V");
    writeDirectors(lcDirector.getDirectors(), "_TNq0_"+toString(dc_volt)+"V");
}

void testTNWithQ0Dynamic(){
    std::cout << "Start synamic TN Calculation with q0 and total time (ms) =" + toString(1000) << std::endl;
    size_t layerNum = 40;
    double dt = 0.01; //ms
    size_t maxIter = 100000; // to 1000ms
    LCD1D::LCParamters lcParam = {4.0, 12, 3.6, 60, 12.0, 6.5, 15.0,  2.0*M_PI/70.0};
    LCD1D::DielecParameters tftpiParam = {0.1, 3.5};
    LCD1D::DielecParameters cfpiParam = {0.1, 3.5};
    LCD1D::RubbingCondition rubbing = {89.0*M_PI/180.0, 45.0*M_PI/180.0, 89.0*M_PI/180.0, 90.0*M_PI/180.0};
    LCD1D::Epsilon epsilonr(layerNum, lcParam.thick);
    epsilonr.setTFTPI(tftpiParam);
    epsilonr.setCFPI(cfpiParam);
    LCD1D::Potential potential(layerNum+1, epsilonr, lcParam.thick/layerNum);

    //create step voltage profiles
    std::map<double, double> profile; //ms, voltage
    profile[0.0] = 2.5;
    profile[200.0] = 5;
    profile[400.0] = -2.5;
    profile[600.0] = -5;
    std::shared_ptr<LCD::WaveformBase> waveform;
    waveform.reset(new LCD::StepWaveform(profile, 800.0));
    potential.createDynamicUpdatePolicy(waveform);

    LCD1D::LCDirector lcDirector(layerNum, lcParam, rubbing, epsilonr);
    lcDirector.createVectorFormUpdater(potential, dt);
    double residual = std::numeric_limits<double>::max();
    size_t iterNum = 0;
    //create dump time steps
    std::set<size_t> dumpSteps;
    dumpSteps.insert(19999); //just before 200ms
    dumpSteps.insert(39999); //just before 400ms
    dumpSteps.insert(59999); //just before 600ms
    dumpSteps.insert(79999); //just before 800ms
    dumpSteps.insert(99999); //just before 1000ms
    while(iterNum <= maxIter){
        potential.update(dt*iterNum);
        lcDirector.update();
        ++iterNum;
        if (dumpSteps.find(iterNum) != dumpSteps.end()){
            writePotentials(potential.getPotentials(), "_TNq0_"+toString(dt*iterNum)+"ms");
            writeDirectors(lcDirector.getDirectors(), "_TNq0_"+toString(dt*iterNum)+"ms");
        }
    };
}


int main(int argc, const char *argv[])
{
    testTNCalculation(2.0);
    testTNCalculation(5.0);
    testTNWithQ0Calculation(2.5);
    testTNWithQ0Calculation(5.0);
    testTNWithQ0Dynamic();
    return 0;
}
