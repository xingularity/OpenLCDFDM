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

#include "LCD1D_LCD1DMain.hpp"
#include <limits>
#include <fstream>
#include <cmath>

//LCDOptics::LIGHTSPECTRUMDATA is std::map<double, double>
LCDOptics::LIGHTSPECTRUMDATA ReadLightSourceSpectrum(std::string _filename){
    LCDOptics::LIGHTSPECTRUMDATA lightSrc;
    std::fstream file;
    file.open(_filename.c_str(),std::fstream::in);
    if (!file){
        std::cout << "Open file failed in ReadLightSourceSpectrum" << std::endl;
        assert(false);
    }
    std::string line;
    while(std::getline(file, line)){
        double lambda, power;
        std::string temp;
        std::stringstream ss(line);
        std::getline(ss, temp, ',');
        lambda = strToNumeric<double>(temp);
        std::getline(ss, temp, ',');
        power = strToNumeric<double>(temp);
        lightSrc[lambda] = power;
    }
    return lightSrc;
}

//LCDOptics::NKData is std::map<double, std::complex<double> >
LCDOptics::NKData ReadIsotropicNKData(std::string _filename){
    LCDOptics::NKData nkSpectrum;
    std::fstream file;
    file.open(_filename.c_str(),std::fstream::in);
    if (!file){
        std::cout << "Open file failed in ReadIsotropicNKData" << std::endl;
        assert(false);
    }
    std::string line;
    while(std::getline(file, line)){
        double lambda, n, k;
        std::string temp;
        std::stringstream ss(line);
        std::getline(ss, temp, ',');
        lambda = strToNumeric<double>(temp);
        std::getline(ss, temp, ',');
        n = strToNumeric<double>(temp);
        std::getline(ss, temp, ',');
        k = strToNumeric<double>(temp);
        nkSpectrum[lambda] = LCD::COMPD(n,-1.0*std::abs(k));
    }
    return nkSpectrum;
}

std::map<double, std::vector<std::complex<double> > > ReadUniaxialNKData(std::string _filename){
    std::map<double, std::vector<std::complex<double> > > nkSpectrum;
    std::fstream file;
    file.open(_filename.c_str(),std::fstream::in);
    if (!file){
        std::cout << "Open file failed in ReadIsotropicNKData" << std::endl;
        assert(false);
    }
    std::string line;
    while(std::getline(file, line)){
        double lambda, no, ko, ne, ke;
        std::vector<std::complex<double> > none(2);
        std::string temp;
        std::stringstream ss(line);
        std::getline(ss, temp, ',');
        lambda = strToNumeric<double>(temp);
        std::getline(ss, temp, ',');
        no = strToNumeric<double>(temp);
        std::getline(ss, temp, ',');
        ko = strToNumeric<double>(temp);
        std::getline(ss, temp, ',');
        ne = strToNumeric<double>(temp);
        std::getline(ss, temp, ',');
        ke = strToNumeric<double>(temp);
        none[0] = LCD::COMPD(no,-1.0*std::abs(ko));
        none[1] = LCD::COMPD(ne,-1.0*std::abs(ke));
        nkSpectrum[lambda] = none;
    }
    return nkSpectrum;
}

void writeTransmissions(LCD1D::TRANSRESULT& data, std::vector<std::vector<std::pair<double, double> > > inAngles,
    std::string output_prefix=""){
    std::ofstream output((output_prefix + ".csv").c_str(), std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < data.size(); ++i)
        for(int j = 0; j < data[i].size(); ++j)
            output << (inAngles[i][j].first)*180.0/M_PI << ", " << (inAngles[i][j].second)*180.0/M_PI << ", " << data[i][j] << std::endl;
    output.close();
}

void writeDirectors(LCD::DOUBLEARRAY2D data, std::string output_prefix=""){
    std::ofstream output((output_prefix + "_Directors.csv").c_str(), std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < data.size(); ++i)
        output << data[i][0] << ", " << data[i][1] << ", " << data[i][2] << std::endl;
    output.close();
}

void writeNormalTransmission(LCD::DOUBLEARRAY1D data, std::string output_prefix=""){
    std::ofstream output((output_prefix + "_NormalTrans.csv").c_str(), std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < data.size(); ++i)
        output << data[i] << std::endl;
    output.close();
}

void testCrossPolarizerNoLC(){
    LCDOptics::LIGHTSPECTRUMDATA lightSrcSpectrum = ReadLightSourceSpectrum("TestLightSrc.csv");
    std::map<double, std::vector<std::complex<double> > > polarizerSpectrum = ReadUniaxialNKData("TestPolarizerSpectrum.csv");
    LCD1D::LCD1DStaticMain lcd1dstaticmain; //use default constructor, no LC calculation.
    lcd1dstaticmain.setOMPThreadNum(1);
    //optical axis
    LCD::DOUBLEARRAY2D axis(1);
    axis[0] = LCD::DOUBLEARRAY1D(2);
    axis[0][0] = 90.0*M_PI/180.0;
    axis[0][1] = 135.0*M_PI/180.0;
    lcd1dstaticmain.addOpticalPolarizer(20.0, polarizerSpectrum, axis);
    axis[0][1] = 45.0*M_PI/180.0;
    lcd1dstaticmain.addOpticalPolarizer(20.0, polarizerSpectrum, axis);
    lcd1dstaticmain.setOpticalIncidentAngles(1,1);
    lcd1dstaticmain.setOpticalWavelength(0.55);
    lcd1dstaticmain.createExtendedJones();
    lcd1dstaticmain.calculate();
    std::vector<LCD1D::TRANSRESULT> trans = lcd1dstaticmain.getTransmissions();
    std::vector<std::vector<std::pair<double, double> > > inAngles = lcd1dstaticmain.getIncidentAngles();
    writeTransmissions(trans[0], inAngles, "CrossPolarizerNoLC_0.55um");
    lcd1dstaticmain.setOpticalWavelength(0.38, 0.78, 0.01);
    lcd1dstaticmain.setOpticalSourceSpectrum(lightSrcSpectrum);
    lcd1dstaticmain.createExtendedJones();
    lcd1dstaticmain.calculate();
    trans = lcd1dstaticmain.getTransmissions();
    writeTransmissions(trans[0], inAngles, "CrossPolarizerNoLC_MultiWavelength");
    lcd1dstaticmain.useOptical2X2Lambertian();
    lcd1dstaticmain.createExtendedJones();
    lcd1dstaticmain.calculate();
    trans = lcd1dstaticmain.getTransmissions();
    writeTransmissions(trans[0], inAngles, "CrossPolarizerNoLC_MultiWavelength_Lambertian");
}

void testTN(){
    LCDOptics::LIGHTSPECTRUMDATA lightSrcSpectrum = ReadLightSourceSpectrum("TestLightSrc.csv");
    std::map<double, std::vector<std::complex<double> > > polarizerSpectrum = ReadUniaxialNKData("TestPolarizerSpectrum.csv");
    std::map<double, std::vector<std::complex<double> > > LCSpectrum = ReadUniaxialNKData("TestLCSpectrum.csv");
    LCD1D::LCParamters LCParam = {4.0, 12, 3.6, 60, 12.0, 6.5, 15.0,  2.0*M_PI/70.0};
    LCD1D::DielecParameters tftpiParam = {0.1, 3.6};
    LCD1D::DielecParameters cfpiParam = {0.1, 3.6};
    LCD1D::RubbingCondition rubbing = {89.0*M_PI/180.0, 45.0*M_PI/180.0, 89.0*M_PI/180.0, 90.0*M_PI/180.0};
    LCD1D::LCD1DStaticMain lcd1dstaticmain(40, 0.01, LCParam, rubbing, 0.0, 7.0, 0.1, 100000000, 1.0e-8);
    lcd1dstaticmain.setTFTPI(tftpiParam);
    lcd1dstaticmain.setCFPI(cfpiParam);

    LCD::DOUBLEARRAY2D axis(1);
    axis[0] = LCD::DOUBLEARRAY1D(2);
    axis[0][0] = 90.0*M_PI/180.0;
    axis[0][1] = 135.0*M_PI/180.0;
    lcd1dstaticmain.addOpticalPolarizer(20.0, polarizerSpectrum, axis);
    lcd1dstaticmain.addOpticalLC(4.8, LCSpectrum);
    axis[0][1] = 45.0*M_PI/180.0;
    lcd1dstaticmain.addOpticalPolarizer(20.0, polarizerSpectrum, axis);
    lcd1dstaticmain.setOpticalIncidentAngles(1,1);
    lcd1dstaticmain.setOpticalWavelength(0.38, 0.78, 0.01);
    lcd1dstaticmain.setOpticalSourceSpectrum(lightSrcSpectrum);
    lcd1dstaticmain.useOptical2X2Lambertian(true);
    lcd1dstaticmain.setOMPThreadNum(8);
    lcd1dstaticmain.createExtendedJones();
    lcd1dstaticmain.calculate();
    std::vector<LCD1D::TRANSRESULT> trans = lcd1dstaticmain.getTransmissions();
    std::vector<std::vector<std::pair<double, double> > > inAngles = lcd1dstaticmain.getIncidentAngles();
    writeTransmissions(trans[50], inAngles, "TestTN_5V_Multi_Lambertian");
    writeTransmissions(trans[20], inAngles, "TestTN_2V_Multi_Lambertian");
    std::vector<LCD::DOUBLEARRAY2D> directors = lcd1dstaticmain.getLCDirResults();
    writeDirectors(directors[20], "TestTN_2V_Multi_Lambertian");
    writeDirectors(directors[50], "TestTN_5V_Multi_Lambertian");
    LCD::DOUBLEARRAY1D normalTrans = lcd1dstaticmain.getNormalTransmissions();
    writeNormalTransmission(normalTrans, "TestTN_Normal_Multi_Lambertian");
}

void testTNDynamic(){
    LCDOptics::LIGHTSPECTRUMDATA lightSrcSpectrum = ReadLightSourceSpectrum("TestLightSrc.csv");
    std::map<double, std::vector<std::complex<double> > > polarizerSpectrum = ReadUniaxialNKData("TestPolarizerSpectrum.csv");
    std::map<double, std::vector<std::complex<double> > > LCSpectrum = ReadUniaxialNKData("TestLCSpectrum.csv");
    LCD1D::LCParamters LCParam = {4.0, 12, 3.6, 60, 12.0, 6.5, 15.0,  2.0*M_PI/70.0};
    LCD1D::DielecParameters tftpiParam = {0.1, 3.6};
    LCD1D::DielecParameters cfpiParam = {0.1, 3.6};
    LCD1D::RubbingCondition rubbing = {89.0*M_PI/180.0, 45.0*M_PI/180.0, 89.0*M_PI/180.0, 90.0*M_PI/180.0};
    LCD1D::LCD1DDynamicMain lcd1ddynamicmain(40, 0.01, LCParam, rubbing, 1000.0);
    lcd1ddynamicmain.setTFTPI(tftpiParam);
    lcd1ddynamicmain.setCFPI(cfpiParam);

    //setup voltage waveform
    std::map<double, double> stepVoltProfile;
    stepVoltProfile[0] = 2.0;
    stepVoltProfile[200] = 5.0;
    stepVoltProfile[400] = -2.0;
    stepVoltProfile[600] = -5.0;
    lcd1ddynamicmain.setStepWaveform(stepVoltProfile, 800.0);

    //setup dump time
    LCD::DOUBLEARRAY1D timeToRecord;
    timeToRecord.push_back(0);
    timeToRecord.push_back(199);
    timeToRecord.push_back(399);
    timeToRecord.push_back(599);
    timeToRecord.push_back(799);
    timeToRecord.push_back(999);
    lcd1ddynamicmain.setRecordTime(timeToRecord);

    //setup optical parameters
    LCD::DOUBLEARRAY2D axis(1);
    axis[0] = LCD::DOUBLEARRAY1D(2);
    axis[0][0] = 90.0*M_PI/180.0;
    axis[0][1] = 135.0*M_PI/180.0;
    lcd1ddynamicmain.addOpticalPolarizer(20.0, polarizerSpectrum, axis);
    lcd1ddynamicmain.addOpticalLC(4.8, LCSpectrum);
    axis[0][1] = 45.0*M_PI/180.0;
    lcd1ddynamicmain.addOpticalPolarizer(20.0, polarizerSpectrum, axis);
    lcd1ddynamicmain.setOpticalIncidentAngles(1,1);
    lcd1ddynamicmain.setOpticalWavelength(0.38, 0.78, 0.01);
    lcd1ddynamicmain.setOpticalSourceSpectrum(lightSrcSpectrum);
    lcd1ddynamicmain.useOptical2X2Lambertian(true);
    lcd1ddynamicmain.setOMPThreadNum(8);
    lcd1ddynamicmain.createExtendedJones();
    lcd1ddynamicmain.calculate();

    std::vector<LCD1D::TRANSRESULT> trans = lcd1ddynamicmain.getTransmissions();
    std::vector<std::vector<std::pair<double, double> > > inAngles = lcd1ddynamicmain.getIncidentAngles();
    std::vector<LCD::DOUBLEARRAY2D> directors = lcd1ddynamicmain.getLCDirResults();
    LCD::DOUBLEARRAY1D recordTime = lcd1ddynamicmain.getRecordTime();

    for (int i = 0; i < recordTime.size(); ++i){
        writeTransmissions(trans[i], inAngles, "TestDyanamicTN_" + toString(recordTime[i]) + "ms_Multi_Lambertian");
        writeDirectors(directors[i], "TestDyanamicTN_" + toString(recordTime[i]) + "ms");
    }
    LCD::DOUBLEARRAY1D normalTrans = lcd1ddynamicmain.getNormalTransmissions();
    writeNormalTransmission(normalTrans, "TestTNDynamic_Normal_Multi_Lambertian");
}

int main(int argc, const char *argv[])
{
    testCrossPolarizerNoLC();
    testTN();
    testTNDynamic();
    return 0;
}
