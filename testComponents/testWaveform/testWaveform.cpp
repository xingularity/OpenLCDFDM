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
 *   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
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

#include "LCD_TimeWaveform.hpp"

using namespace LCD;

void testDC(){
    std::cout << "test DCWaveform" << std::endl;
    DCWaveform dc(5);
    std::cout << "10 ms: " << dc(10.0) << std::endl;
    DCWaveform dc2(-5.2, 10.2, 3.0);
    std::cout << "10 ms: " << dc2(10.0) << ", 12ms:" << dc2(12)<< std::endl;
}

void testStep(){
    std::cout << "test StepWaveform" << std::endl;
    std::map<double, double> input;
    input[12.3] = 5.0;
    input[15] = -4.5;
    input[19.3] = 0.3;
    input[21.5] = -0.5;
    double wavePeriod = 25.0;//ms
    StepWaveform step(input, wavePeriod);
    std::cout << "print profile: " << std::endl;
    for (auto i : step.getStepProfile()){
        std::cout << i[0] << ": " << i[1] << std::endl;
    }
    std::cout << "No ampShift and no pre apply" << std::endl;
    std::cout << "0 ms: " << step(0.0) << ", 2.8ms:" << step(2.8) << ", 3.3ms: " << step(3.3) << ", 7.2ms: " << step(7.2)
    << ", 9.7ms: " << step(9.7) << ", 12.7ms: " << step(12.7) << ", 24.9ms: " << step(24.9) << ", 27.8ms: "<< step(27.8) << std::endl;
    std::cout << "ampShift=2.0 and preApplyTime=10, preApplyVal=3.0" << std::endl;
    StepWaveform step2(input, wavePeriod, 2.0, 10.0, 3.0);
    std::cout << "0 ms: " << step2(0.0) << ", 5 ms: " << step2(5)  << ", 12.6ms:" << step2(12.6) << ", 13.3ms: " << step2(13.3) << ", 17.2ms: " << step2(17.2)
    << ", 19.7ms: " << step2(19.7) << ", 22.7ms: " << step2(22.7) << ", 34.9ms: " << step2(34.9) << ", 37.6ms: "<< step2(37.6) << ", 60.5 ms: " << step2(60.5) << std::endl;
}

int main(int argc, const char *argv[])
{
    testDC();
    testStep();
    return 0;
}
