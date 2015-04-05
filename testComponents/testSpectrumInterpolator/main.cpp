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

#include <iostream>
#include <list>
#include "../../include/LCD_OpticsBaseDef.hpp"

using namespace std;
using namespace LCDOptics;

int main(int argc, const char *argv[])
{
    //give spectru data
    LIGHTSPECTRUMDATA lightSrc;
    lightSrc[0.38] = 1.2;
    lightSrc[0.63] = -3.4;
    lightSrc[0.78] = 3.6;
    NKData nkData;
    nkData[0.38] = COMPD(1.2, -3.2);
    nkData[0.63] = COMPD(0.9, 3.5);
    nkData[0.78] = COMPD(-2.1, 8.4);
    NKoNKeData nkonkeData;
    nkonkeData[0.38] = std::make_tuple(COMPD(1.2, -3.6), COMPD(1.2, -2.3));
    nkonkeData[0.63] = std::make_tuple(COMPD(2.4, 3.6), COMPD(3.6, -4.2));
    nkonkeData[0.78] = std::make_tuple(COMPD(3.6, -3.6), COMPD(2.4, 3.2));
    NK123Data nk123Data;
    nk123Data[0.38] = std::make_tuple(COMPD(1.2, -1.8), COMPD(1.5, 0.0), COMPD(2.0, -4.0));
    nk123Data[0.63] = std::make_tuple(COMPD(2.4, 3.6), COMPD(3.0, 3.2), COMPD(3.6, -1.0));
    nk123Data[0.78] = std::make_tuple(COMPD(3.6, 2.4), COMPD(4.5, -3.6), COMPD(1.8, -3.0));
    DOUBLEARRAY1D lambdas;
    lambdas.push_back(0.35);
    lambdas.push_back(0.55);
    lambdas.push_back(0.63);
    lambdas.push_back(0.75);
    lambdas.push_back(0.81);

    //This line is to test static assert.
    //SpectrumInterpolator<std::map<COMPD, double> > testStaticAssert(lambdas);

    SpectrumInterpolator<LIGHTSPECTRUMDATA> lightInterp(lambdas);
    std::vector<LIGHTSPECTRUMDATA::value_type::second_type> lightOut;
    lightInterp.interpolate(lightSrc, lightOut);
    cout << "=====type: LIGHTSPECTRUMDATA=====" <<endl;
    for (int i =0; i < lightOut.size(); ++i){
        cout << lambdas[i] << ":" << lightOut[i] <<endl;
    }
    cout << endl;
    cout << "0.55: " << lightInterp.interpolate(lightSrc, 0.55) << endl;
    cout << "0.6: " << lightInterp.interpolate(lightSrc, 0.6) << endl;
    cout << "0.9: " << lightInterp.interpolate(lightSrc, 0.9) << endl;

    SpectrumInterpolator<NKData> isoInterp(lambdas);
    std::vector<NKData::value_type::second_type> isoOut;
    isoInterp.interpolate(nkData, isoOut);
    cout << "=====type: NKData=====" <<endl;
    for (int i =0; i < isoOut.size(); ++i){
        cout << lambdas[i] << ":" << isoOut[i] <<endl;
    }
    cout << endl;
    cout << "0.55: " << isoInterp.interpolate(nkData, 0.55) << endl;
    cout << "0.6: " << isoInterp.interpolate(nkData, 0.6) << endl;
    cout << "0.9: " << isoInterp.interpolate(nkData, 0.9) << endl;

    SpectrumInterpolator<NKoNKeData> uniInterp(lambdas);
    std::vector<NKoNKeData::value_type::second_type> uniOut;
    uniInterp.interpolate(nkonkeData, uniOut);
    cout << "=====type: NKoNKeData=====" <<endl;
    for (int i =0; i < uniOut.size(); ++i){
        cout << lambdas[i] << ":" << get<0>(uniOut[i]) << ", " << get<1>(uniOut[i]) <<endl;
    }
    cout << endl;
    cout << "0.55: " << get<0>(uniInterp.interpolate(nkonkeData, 0.55)) << ", " << get<1>(uniInterp.interpolate(nkonkeData, 0.55)) << endl;
    cout << "0.6: " << get<0>(uniInterp.interpolate(nkonkeData, 0.6)) << ", " << get<1>(uniInterp.interpolate(nkonkeData, 0.6)) << endl;
    cout << "0.9: " << get<0>(uniInterp.interpolate(nkonkeData, 0.9)) << ", " << get<1>(uniInterp.interpolate(nkonkeData, 0.9)) << endl;

    SpectrumInterpolator<NK123Data> biInterp(lambdas);
    std::vector<NK123Data::value_type::second_type> biOut;
    biInterp.interpolate(nk123Data, biOut);
    cout << "=====type: NK123Data=====" <<endl;
    for (int i =0; i < biOut.size(); ++i){
        cout << lambdas[i] << ":" << get<0>(biOut[i]) << ", " << get<1>(biOut[i]) << ", " << get<2>(biOut[i]) <<endl;
    }
    cout << endl;
    cout << "0.55: " << get<0>(biInterp.interpolate(nk123Data, 0.55)) << ", " << get<1>(biInterp.interpolate(nk123Data, 0.55)) << ", " << get<2>(biInterp.interpolate(nk123Data, 0.55)) << endl;
    cout << "0.6: " << get<0>(biInterp.interpolate(nk123Data, 0.6)) << ", " << get<1>(biInterp.interpolate(nk123Data, 0.6)) << ", " << get<2>(biInterp.interpolate(nk123Data, 0.6)) << endl;
    cout << "0.9: " << get<0>(biInterp.interpolate(nk123Data, 0.9)) << ", " << get<1>(biInterp.interpolate(nk123Data, 0.9)) << ", " << get<2>(biInterp.interpolate(nk123Data, 0.9)) << endl;

    return 0;
}
