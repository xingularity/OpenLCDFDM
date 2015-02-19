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
        nkSpectrum[lambda] = COMPD(n,k);
    }
    return nkSpectrum;
}

void testSingleGlass(){
    std::cout << "Start to calculate single glass case..." << std::endl;
    NKData nkSpectrum = ReadIsotropicNKData("GlassNKSpectrum.csv");
    LIGHTSPECTRUMDATA lightSrc = ReadLightSourceSpectrum("D65.csv");
    /*
    for (auto i : nkSpectrum){
        std:: cout << i.first << ", " << i.second << std::endl;
    }
    */
    /*
    for (auto i : lightSrc){
        std:: cout << i.first << ", " << i.second << std::endl;
    }
    */
    

}

int main(int argc, char const *argv[]) {
    testSingleGlass();
    return 0;
}
