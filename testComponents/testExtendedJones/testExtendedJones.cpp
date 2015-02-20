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

#include <iostream>
#include <fstream>
#include <iomanip>
#include "../../include/LCD1D_ExtendedJones.hpp"
#include "../../include/LCD_UsefulFuncs.hpp"

using namespace std;
using namespace LCDOptics;
using namespace LCD1D;
using namespace LCD;

IAngles createAngles(){
    IAngles inAngles;

}

LIGHTSPECTRUMDATA ReadLightSourceSpectrum(std::string _filename){
    LIGHTSPECTRUMDATA lightSrc;
    fstream file;
    file.open(_filename.c_str(),fstream::in);
    if (!file){
        std::cout << "Open file failed in ReadLightSourceSpectrum" << std::endl;
        assert(false);
    }
    string line;
    while(getline(file, line)){
        double lambda, power;
        std::string temp;
        stringstream ss(line);
        getline(ss, temp, ',');
        lambda = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        power = strToNumeric<double>(temp);
        lightSrc[lambda] = power;
    }
    return lightSrc;
}

NKData ReadIsotropicNKData(std::string _filename){
    NKData nkSpectrum;
    fstream file;
    file.open(_filename.c_str(),fstream::in);
    if (!file){
        std::cout << "Open file failed in ReadIsotropicNKData" << std::endl;
        assert(false);
    }
    string line;
    while(getline(file, line)){
        double lambda, n, k;
        std::string temp;
        stringstream ss(line);
        getline(ss, temp, ',');
        lambda = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        n = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        k = strToNumeric<double>(temp);
        nkSpectrum[lambda] = COMPD(n,-1.0*std::abs(k));
    }
    return nkSpectrum;
}

//create all viweing angles
IAngles createInAngles(){
    IAngles answer;
    answer = std::vector< std::vector<Angle> >(81, std::vector<Angle>(361));
    for (int i = 0; i <= 80; ++i)
        for(int j = 0; j <= 360; ++j)
            answer[i][j] = makeAngle2(i,j);
    return answer;
}

void testSingleGlassMultiWavelength(){
    std::cout << "Start to calculate single glass case..." << std::endl;
    NKData nkSpectrum = ReadIsotropicNKData("GlassNKSpectrum.csv");
    LIGHTSPECTRUMDATA lightSrc = ReadLightSourceSpectrum("D65.csv");
    MATERIALLAYERS2X2CONT materials;
    materials.push_back(Optical2x2IsoPtr(new Optical2X2OneLayer<ISOType>(500.0, nkSpectrum, OPT_GLASS)));
    IAngles inAngles = createInAngles();
    ExtendedJones extj(materials, inAngles, 0.38, 0.78, 0.01, lightSrc, false);
    extj.calculateExtendedJones();
    TRANSRESULT answer = extj.getTransmissions();
    //output to file
    ofstream output("GlassOnly_MultiWaveLength_D65.csv", std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < answer.size(); ++i)
        for(int j = 0; j < answer[i].size(); ++j)
            output << std::get<0>(inAngles[i][j])*180.0/M_PI << ", " << std::get<1>(inAngles[i][j])*180.0/M_PI << "," << answer[i][j] << std::endl;
    output.close();
}

void testSingleGlassSingleWavelength(){
    std::cout << "Start to calculate single glass case..." << std::endl;
    NKData nkSpectrum = ReadIsotropicNKData("GlassNKSpectrum.csv");
    LIGHTSPECTRUMDATA lightSrc = ReadLightSourceSpectrum("D65.csv");
    MATERIALLAYERS2X2CONT materials;
    materials.push_back(Optical2x2IsoPtr(new Optical2X2OneLayer<ISOType>(500.0, nkSpectrum, OPT_GLASS)));
    IAngles inAngles = createInAngles();
    ExtendedJones extj(materials, inAngles, 0.55, lightSrc, false);
    extj.calculateExtendedJones();
    TRANSRESULT answer = extj.getTransmissions();
    //output to file
    ofstream output("GlassOnly_SingleWaveLength.csv", std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < answer.size(); ++i)
        for(int j = 0; j < answer[i].size(); ++j)
            output << std::get<0>(inAngles[i][j])*180.0/M_PI << ", " << std::get<1>(inAngles[i][j])*180.0/M_PI << "," << answer[i][j] << std::endl;
    output.close();
}

int main(int argc, char const *argv[]) {
    testSingleGlassSingleWavelength();
    testSingleGlassMultiWavelength();
    return 0;
}
