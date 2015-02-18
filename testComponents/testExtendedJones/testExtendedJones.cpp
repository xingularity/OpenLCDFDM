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
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace LCDOptics;
using namespace LCD1D;
using namespace LCD;

NKData ReadIsotropicNKData(std::string _filename){
    NKData nkSpectrum;
    fstream file;
    file.open(_filename.c_str(),fstream::in);
    string line;
    while(getline(file, line)){
        double lambda, n, k;
        stringstream ss;
        ss << line;
        ss >> lambda;
        ss >> n;
        ss >> k;
        nkSpectrum[lambda] = COMPD(n,k);
        std::cout << lambda << ", " << nkSpectrum[lambda] << std::endl;
    }
    return nkSpectrum;
}

void testSingleGlass(){
    std::cout << "Start to calculate single glass case..." << std::endl;
    NKData nkSpectrum = ReadIsotropicNKData("GlassNKSpectrum.csv");
    for (auto i : nkSpectrum){
        std:: cout << i.first << ", " << i.second << std::endl;
    }
}

int main(int argc, char const *argv[]) {
    testSingleGlass();
    return 0;
}
