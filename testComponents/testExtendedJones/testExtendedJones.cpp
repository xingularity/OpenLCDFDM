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
#include <omp.h>
#include <ctime>
#include "../../include/LCD1D_ExtendedJones.hpp"
#include "../../include/LCD_UsefulFuncs.hpp"

using namespace std;
using namespace LCDOptics;
using namespace LCD1D;
using namespace LCD;

int ompNumThreads = 0;

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

NKoNKeData ReadUniaxialNKData(std::string _filename){
    NKoNKeData nkSpectrum;
    fstream file;
    file.open(_filename.c_str(),fstream::in);
    if (!file){
        std::cout << "Open file failed in ReadIsotropicNKData" << std::endl;
        assert(false);
    }
    string line;
    while(getline(file, line)){
        double lambda, no, ko, ne, ke;
        std::string temp;
        stringstream ss(line);
        getline(ss, temp, ',');
        lambda = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        no = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        ko = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        ne = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        ke = strToNumeric<double>(temp);
        nkSpectrum[lambda] = std::make_tuple(COMPD(no,-1.0*std::abs(ko)), COMPD(ne,-1.0*std::abs(ke)));
    }
    return nkSpectrum;
}

void ReadLCDirectors(std::string _filename, DIRVEC& dirVec){
    fstream file;
    file.open(_filename.c_str(),fstream::in);
    std::string line;
    //check how many lines in file
    int numLines = 0;
    while(getline(file, line)){
        numLines++;
    }
    file.close();
    //reopen to read
    file.open(_filename.c_str(),fstream::in);
    DIRVEC input;
    input.resize(numLines);
    for(int i = 0; getline(file, line); ++i){
        double nx, ny, nz;
        std::string temp;
        stringstream ss(line);
        getline(ss, temp, ',');
        input(i)(0) = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        input(i)(1) = strToNumeric<double>(temp);
        getline(ss, temp, ',');
        input(i)(2) = strToNumeric<double>(temp);
    }
    file.close();
    dirVec.resize(input.size()-1);
    for(int i =0; i < dirVec.size(); ++i){
        double theta, phi;
        theta = 0.5*(acos(input(i)(2)) + acos(input(i+1)(2)));
        phi = 0.5*(atan2(input(i)(1), input(i)(0)) + atan2(input(i+1)(1), input(i+1)(0)));
        dirVec(i)(0) = sin(theta)*cos(phi);
        dirVec(i)(1) = sin(theta)*sin(phi);
        dirVec(i)(2) = cos(theta);
    }
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

void calculateAndOutput(std::string output_prefix, MATERIALLAYERS2X2CONT materials, LIGHTSPECTRUMDATA lightSrc){
    IAngles inAngles = createInAngles();
    //start to calculate single wavelength case
    ExtendedJones extj(materials, inAngles, 0.55, lightSrc, false, false);
    extj.calculateExtendedJones();
    TRANSRESULT answer = extj.getTransmissions();
    //output to file
    ofstream output((output_prefix + "_SingleWaveLength_0.55um.csv").c_str(), std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < answer.size(); ++i)
        for(int j = 0; j < answer[i].size(); ++j)
            output << std::get<0>(inAngles[i][j])*180.0/M_PI << ", " << std::get<1>(inAngles[i][j])*180.0/M_PI << ", " << answer[i][j] << std::endl;
    output.close();

    //start to calculate multiple wavelength case
    ExtendedJones extj2(materials, inAngles, 0.38, 0.78, 0.01, lightSrc, false, false);
    extj2.calculateExtendedJones();
    answer = extj2.getTransmissions();
    //output to file
    output.open((output_prefix + "_MultiWaveLength_TestLightSrc.csv").c_str(), std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < answer.size(); ++i)
        for(int j = 0; j < answer[i].size(); ++j)
            output << std::get<0>(inAngles[i][j])*180.0/M_PI << ", " << std::get<1>(inAngles[i][j])*180.0/M_PI << ", " << answer[i][j] << std::endl;
    output.close();

    //start to calculate multiple wavelength case
    ExtendedJones extj3(materials, inAngles, 0.38, 0.78, 0.01, LIGHTSPECTRUMDATA(), false, false);
    extj3.calculateExtendedJones();
    answer = extj3.getTransmissions();
    //output to file
    output.open((output_prefix + "_MultiWaveLength_EqWhite.csv").c_str(), std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < answer.size(); ++i)
        for(int j = 0; j < answer[i].size(); ++j)
            output << std::get<0>(inAngles[i][j])*180.0/M_PI << ", " << std::get<1>(inAngles[i][j])*180.0/M_PI << ", " << answer[i][j] << std::endl;
    output.close();

    //start to calculate multiple wavelength case with Lambertian light source
    ExtendedJones extj4(materials, inAngles, 0.38, 0.78, 0.01, lightSrc, true, false);
    extj4.calculateExtendedJones();
    answer = extj4.getTransmissions();
    //output to file
    output.open((output_prefix + "_MultiWaveLength_TestLightSrc_Lambertian.csv").c_str(), std::fstream::out|std::fstream::trunc);
    output << std::setprecision(15);
    for (int i =0; i < answer.size(); ++i)
        for(int j = 0; j < answer[i].size(); ++j)
            output << std::get<0>(inAngles[i][j])*180.0/M_PI << ", " << std::get<1>(inAngles[i][j])*180.0/M_PI << ", " << answer[i][j] << std::endl;
    output.close();
}

void testSingleGlass(){
    std::cout << "Start to calculate single glass case..." << std::endl;
    NKData nkSpectrum = ReadIsotropicNKData("TestGlassNKSpectrum.csv");
    LIGHTSPECTRUMDATA lightSrc = ReadLightSourceSpectrum("TestLightSrc.csv");
    MATERIALLAYERS2X2CONT materials;
    materials.push_back(Optical2x2IsoPtr(new Optical2X2OneLayer<ISOType>(500.0, nkSpectrum, OPT_GLASS)));
    calculateAndOutput("GlassOnly", materials, lightSrc);
}

void testCrossPolarizer(){
    std::cout << "Start to calculate cross polarizer case..." << std::endl;
    LIGHTSPECTRUMDATA lightSrc = ReadLightSourceSpectrum("TestLightSrc.csv");
    NKoNKeData polarizerSpectrum = ReadUniaxialNKData("TestPolarizerSpectrum.csv");
    MATERIALLAYERS2X2CONT materials;
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(20.0, polarizerSpectrum, OPT_POLARIZER)));
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(20.0, polarizerSpectrum, OPT_POLARIZER)));
    //set first polarizer to 45 degree
    double theta_pol = 90.0*M_PI/180.0;
    double phi_pol1 = 135.0*M_PI/180.0;
    double phi_pol2 = 45.0*M_PI/180.0;
    DIRVEC pol1; pol1.resize(1);
    pol1(0)(0) = sin(theta_pol)*cos(phi_pol1);
    pol1(0)(1) = sin(theta_pol)*sin(phi_pol1);
    pol1(0)(2) = cos(theta_pol);
    materials[0]->resetAxes(pol1);
    DIRVEC pol2; pol2.resize(1);
    pol2(0)(0) = sin(theta_pol)*cos(phi_pol2);
    pol2(0)(1) = sin(theta_pol)*sin(phi_pol2);
    pol2(0)(2) = cos(theta_pol);
    materials[1]->resetAxes(pol2);
    calculateAndOutput("CrossPolarizer", materials, lightSrc);
}

void testLC(){
    std::cout << "Start to calculate LC no glass case..." << std::endl;
    DIRVEC dirVec;
    ReadLCDirectors("TN_Director.txt", dirVec);
    LIGHTSPECTRUMDATA lightSrc = ReadLightSourceSpectrum("TestLightSrc.csv");
    NKoNKeData polarizerSpectrum = ReadUniaxialNKData("TestPolarizerSpectrum.csv");
    NKoNKeData LCSpectrum = ReadUniaxialNKData("TestPositiveLC.csv");
    MATERIALLAYERS2X2CONT materials;
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(20.0, polarizerSpectrum, OPT_POLARIZER)));
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(4.8, LCSpectrum, OPT_LCMATERIAL)));
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(20.0, polarizerSpectrum, OPT_POLARIZER)));
    //set first polarizer to 45 degree
    double theta_pol = 90.0*M_PI/180.0;
    double phi_pol1 = 135.0*M_PI/180.0;
    double phi_pol2 = 45.0*M_PI/180.0;
    DIRVEC pol1; pol1.resize(1);
    pol1(0)(0) = sin(theta_pol)*cos(phi_pol1);
    pol1(0)(1) = sin(theta_pol)*sin(phi_pol1);
    pol1(0)(2) = cos(theta_pol);
    materials[0]->resetAxes(pol1);
    DIRVEC pol2; pol2.resize(1);
    pol2(0)(0) = sin(theta_pol)*cos(phi_pol2);
    pol2(0)(1) = sin(theta_pol)*sin(phi_pol2);
    pol2(0)(2) = cos(theta_pol);
    materials[2]->resetAxes(pol2);
    materials[1]->resetAxes(dirVec);
    calculateAndOutput("TN_NoGlass", materials, lightSrc);
}

void testLCWithGlass(){
    std::cout << "Start to calculate LC with glass case..." << std::endl;
    DIRVEC dirVec;
    ReadLCDirectors("TN_Director.txt", dirVec);
    NKData nkSpectrum = ReadIsotropicNKData("TestGlassNKSpectrum.csv");
    LIGHTSPECTRUMDATA lightSrc = ReadLightSourceSpectrum("TestLightSrc.csv");
    NKoNKeData polarizerSpectrum = ReadUniaxialNKData("TestPolarizerSpectrum.csv");
    NKoNKeData LCSpectrum = ReadUniaxialNKData("TestPositiveLC.csv");
    MATERIALLAYERS2X2CONT materials;
    materials.push_back(Optical2x2IsoPtr(new Optical2x2IsoPtr::element_type(500.0, nkSpectrum, OPT_GLASS)));
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(20.0, polarizerSpectrum, OPT_POLARIZER)));
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(4.8, LCSpectrum, OPT_LCMATERIAL)));
    materials.push_back(Optical2x2UnixialPtr(new Optical2x2UnixialPtr::element_type(20.0, polarizerSpectrum, OPT_POLARIZER)));
    materials.push_back(Optical2x2IsoPtr(new Optical2x2IsoPtr::element_type(500.0, nkSpectrum, OPT_GLASS)));
    //set first polarizer to 45 degree
    double theta_pol = 90.0*M_PI/180.0;
    double phi_pol1 = 135.0*M_PI/180.0;
    double phi_pol2 = 45.0*M_PI/180.0;
    DIRVEC pol1; pol1.resize(1);
    pol1(0)(0) = sin(theta_pol)*cos(phi_pol1);
    pol1(0)(1) = sin(theta_pol)*sin(phi_pol1);
    pol1(0)(2) = cos(theta_pol);
    materials[1]->resetAxes(pol1);
    DIRVEC pol2; pol2.resize(1);
    pol2(0)(0) = sin(theta_pol)*cos(phi_pol2);
    pol2(0)(1) = sin(theta_pol)*sin(phi_pol2);
    pol2(0)(2) = cos(theta_pol);
    materials[3]->resetAxes(pol2);
    materials[2]->resetAxes(dirVec);
    calculateAndOutput("TN_WithGlass", materials, lightSrc);
}

int main(int argc, char const *argv[]) {
    time_t start_time, end_time, diff_time;
    std::time(&start_time);
    std::cout << "omp_get_num_procs() = " << omp_get_num_procs() << std::endl;
    if (argc > 1){
        ompNumThreads = atoi(argv[1]);
    }
    omp_set_num_threads(ompNumThreads);
    testSingleGlass();
    testCrossPolarizer();
    testLC();
    testLCWithGlass();
    std::time(&end_time);
    diff_time = std::difftime(end_time, start_time);
    std::cout << "Total time consumed: " << diff_time << "(s)" << std::endl;
    return 0;
}
