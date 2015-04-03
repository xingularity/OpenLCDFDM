#include "LCD1D_FDM1DSolver.hpp"
#include <limits>
#include <fstream>
#include <cmath>
#include "LCD_UsefulFuncs.hpp"

void writeDirectors(LCD1D::DIRVEC dirs, std::string suffix=""){
    std::fstream file;
    file.open(("Directors"+suffix+".txt").c_str(), std::ios::trunc);
    for (int i = 0; i < dirs.extent(0); ++i){
        file << dirs(i)(0) << " " << dirs(i)(1) << " " << dirs(i)(2) << " " << std::acos(dirs(i)(2))
             << " " << std::atan2(dirs(i)(1), dirs(i)(0)) << " " << std::endl;
    }
    file.flush();
    file.close();
}

void writePotentials(LCD::DOUBLEARRAY1D potentials, std::string suffix=""){
    std::fstream file;
    file.open(("Potentials"+suffix+".txt").c_str(), std::ios::trunc);
    for (int i = 0; i < potentials.size(); ++i){
        file << potentials[i] << std::endl;
    }
    file.flush();
    file.close();
}

void testTNCalculation(double dc_volt){
    std::cout << "Start TN Calculation without q0 under dc_volt=" + toString(dc_volt) << std::endl;
    size_t layerNum = 48;
    double dt = 0.01; //ms
    size_t maxIter = 100000000;
    double max_error = 1.0e-8;
    LCD1D::LCParamters lcParam = {4.8, 11.98, 3.36, 48.0, 11.5, 6, 14.3, 0.0};
    LCD1D::DielecParameters tftpiParam = {0.1, 3.2};
    LCD1D::DielecParameters cfpiParam = {0.1, 3.2};
    LCD1D::RubbingCondition rubbing = {89.0*M_PI/180.0, 45.0*M_PI/180.0, 89.0*M_PI/180.0, 90.0*M_PI/180.0};
    LCD1D::Epsilon epsilonr(48, lcParam.thick);
    epsilonr.setTFTPI(tftpiParam);
    epsilonr.setCFPI(cfpiParam);
    LCD1D::Potential potential(layerNum, epsilonr, lcParam.thick/layerNum);
    potential.createStaticUpdatePolicy();
    LCD1D::LCDirector lcDirector(layerNum, lcParam, rubbing, epsilonr);
    lcDirector.createVectorFormUpdater(potential, dt);
    double residual = std::numeric_limits<double>::max();
    while(residual >= max_error){
        potential.update(dc_volt);
        residual = lcDirector.update();
    };
    writePotentials(potential.getPotentials(), "_TNNoq0_"+toString(dc_volt)+"V");
    writeDirectors(lcDirector.getDirectors(), "_TNNoq0_"+toString(dc_volt)+"V");
}

int main(int argc, const char *argv[])
{
    testTNCalculation(2.0);
    return 0;
}
