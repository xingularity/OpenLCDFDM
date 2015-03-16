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
