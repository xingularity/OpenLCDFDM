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
#include <omp.h>
#include <limits>

LCD1DMainBase::LCD1DMainBase(double _lcLayerNum, double _dt, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing){
	dt = _dt;
	epsilonr.reset(new LCD1D::Epsilon(_lcLayerNum, _lcParam.thick);
	potentials.reset(new LCD1D::Potential(_lcLayerNum+1, *epsilonr, _lcParam.thick/_lcLayerNum));
	lcDir.reset(new LCD1D::LCDirector(_lcLayerNum, _lcParam, _rubbing, *epsilonr));
	lcDirector.createVectorFormUpdater(potential, dt);
}

LCD1DMainBase::LCD1DMainBase(){}

void LCD1DMainBase::setTFTPI(LCD1D::DielecParameters _tftpi){
	epsilonr->setTFTPI(_tftpi);
}

void LCD1DMainBase::setCFPI(LCD1D::DielecParameters _cfpi){
	epsilonr->setCFPI(_cfpi);
}

void LCD1DMainBase::setOMPThreadNum(size_t _num){
	int num = _num;
	if (num > omp_get_num_procs())
		num = omp_get_num_procs();
	omp_set_num_threads(num);
}

size_t LCD1DMainBase::addOpticalGlassLayer(double _thick, std::map<double, std::complex<double> > _nkSpectrum, int pos){
	return addOpticalIsotropicLayer(_thick, _nkSpectrum, LCDOptics::OPT_GLASS, pos);
}

size_t LCD1DMainBase::addOpticalIsotropicLayer(double _thick, std::map<double, std::complex<double> > _nkSpectrum,
    LCDOptics::OpticalMaterialClass _class, int pos){
	if (pos < 0) || (pos >= materials.size()){
		materials.push_back(LCDOptics::Optical2x2IsoPtr(
			new LCDOptics::Optical2x2IsoPtr::element_type(thick, _nkSpectrum, _class)));
		return materials.size()-1;
	}
	else{
		materials.insert(materials.begin() + pos, LCDOptics::Optical2x2IsoPtr(
			new LCDOptics::Optical2x2IsoPtr::element_type(thick, _nkSpectrum, _class)));
		return pos;
	}
}

size_t LCD1DMainBase::addOpticalUnaixialLayer(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum,
    LCDOptics::OpticalMaterialClass _class, LCD::DOUBLEARRAY2D _axes, int pos){

	LCDOptics::NKoNKeData nkSpectrum;
	for (auto& data : _nkSpectrum){
        if (data.second.size() != 2){
            std::string msg = "Size of the nkSpectrum data for uniaxial material at wavelength ";
            msg += toString(data.first);
    l    	msg += " is not 2.";
            throw std::runtime_error(msg.c_str());
        }
        //make sure the k value is negative
		COMPD nko(data.second[0].real(), -1.0*std::abs(data.second[0].imag()));
		COMPD nke(data.second[1].real(), -1.0*std::abs(data.second[1].imag()));
		nkSpectrum[data.first] = std::make_tuple(nko, nke);
	}

	DIRVEC axes;
	axes.resize(axes.size());
	for (int i = 0; i < axes.extent(0); ++i){
		if (_axes[i].size() != 2)
			throw std::runtime_error("incorrect number of angles of optical axis in setting uniaxial layer");
		axes(i)(0) = sin(_axes[i][0])*cos(_axes[i][1]);
		axes(i)(1) = sin(_axes[i][0])*sin(_axes[i][1]);
		axes(i)(2) = cos(_axes[i][0]);
	}

	if (pos < 0) || (pos >= materials.size()){
		materials.push_back(LCDOptics::Optical2x2UnixialPtr(
			new LCDOptics::Optical2x2UnixialPtr::element_type(thick, nkSpectrum, _class)));
		(materials.back())->resetAxes(axes);
		return materials.size() - 1;
	}
	else{
		materials.insert(materials.begin() + pos, LCDOptics::Optical2x2UnixialPtr(
			new LCDOptics::Optical2x2UnixialPtr::element_type(thick, nkSpectrum, _class)));
		(materials[pos])->resetAxes(axes);
		return pos;
	}
}

size_t LCD1DMainBase::addOpticalPolarizer(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum, LCD::DOUBLEARRAY2D _axes, int pos){
    return addOpticalUnaixialLayer(_thick, _nkSpectrum, LCDOptics::OPT_POLARIZER, _axes, pos);
}

size_t LCD1DMainBase::addOpticalLC(double _thick, std::map<double, std::vector<std::complex<double> > > _nkSpectrum, int pos){
    return addOpticalUnaixialLayer(_thick, _nkSpectrum, LCDOptics::OPT_LCMATERIAL, LCD::DOUBLEARRAY2D(0), pos);
}

void LCD1DMainBase::removeOpticalLayer(size_t _index){
	if (_index > materials.size() -1){
		std::string msg = "No such optical layer: " + toString(_index);
		throw std::runtime_error(msg.c_str());
	}
	materials.pop(_index);
}

void setOpticalIncidentAngles(){
	inAngles.clear();
	inAngles = std::vector< std::vector<Angle> >(1, std::vector<Angle>(1));
	inAngles[0][0] = LCDOptics::makeAngle2(0.0, 0.0);
}

void setOpticalIncidentAngles(unsigned double _thetaInterval, unsigned double _phiInterval){
	std::vector<double> thetas(1, 0.0);
    std::vector<double> phis(1, 0.0);

	double acc_error = std::numeric_limits<double>::epsilon()*360.0;

    while(thetas.back()+_thetaInterval <= (360+acc_error)){
        thetas.push_back(thetas.back() + _thetaInterval);
    };

    while(phis.back()+_phiInterval <= (80+acc_error)){
        phis.push_back(phis.back() + _phiInterval);
    };

    inAngles = std::vector< std::vector<Angle> >(thetas.size(), std::vector<Angle>(phis.size()));
    for (int i = 0; i <= thetas.size(); ++i)
        for(int j = 0; j <= phis.size(); ++j)
            inAngles[i][j] = LCDOptics::makeAngle2(thetas[i],phis[j]);

}

void LCD1DMainBase::setOpticalIncidentAngles(std::vector<std::pair<double, double> > _angles){
    inAngles.clear();
    //IAngles is a 2D std::vector structure, I will put all manually input angles in the second dimension.
    inAngles.push_back(1, std::vector<Angle>(0));
    for (auto& i : _angles)
        inAngles[0].push_back(LCDOptics::makeAngle2(i.first, i.second));
}

void LCD1DMainBase::setOpticalWavelength(unsigned double _lambda_start, unsigned double _lambda_end, unsigned double _lambda_step){
	if (_lambda_end < _lambda_start){
		double temp = _lambda_end;
		_lambda_end = _lambda_start;
		_lambda_start = temp;
	}

	if (_lambda_step > (_lambda_end - _lambda_start)){
		_lambda_step = (_lambda_end - _lambda_start)*0.1;
	}

    multiWavelengthLambdas = std::make_tuple(_lambda_start, _lambda_end, _lambda_step);
}

void LCD1DMainBase::setOpticalWavelength(unsigned double _lambda){
    multiWavelengthLambdas = std::make_tuple(_lambda, _lambda, 0);
}

void LCD1DMainBase::setOpticalSourceSpectrum(LCDOptics::LIGHTSPECTRUMDATA _input){
    lightSrcSpectrum = _input;
}

void LCD1DMainBase::enableOptical2X2Calculation(bool _ifDo){
	ifCalculate2X2Optics = _ifDo;
	if (ifCalculate2X2Optics)
		createExtendedJones();
}

void LCD1DMainBase::useOptical2X2Lambertian(bool _if=true){
	ifUseLambertian = _if;
}

void LCD1DMainBase::resetLCParam(const LCD1D:LCParamters _param){
    if (!lcDir){
        std::cout << "No LC calculation, change of LC parameters doesn't happen." << std::endl;
        return
    }

	lcDir->resetConditions(_param);

	if (extJonesMain){
		extJonesMain->resetLCThickness(_param.thick);
	}
}
void LCD1DMainBase::resetLCRubbing(const LCD1D:RubbingCondition _rubbing){
    if (!lcDir){
        std::cout << "No LC calculation, change of LC rubbing condition doesn't happen." << std::endl;
        return
    }
    lcDir->resetConditions(_rubbing);
}

void LCD1DMainBase::createExtendedJones(){
	if (!ifCalculate2X2Optics) return;
	double lambda_start, lambda_end, lambda_step;
	std::tie(lambda_start, lambda_end, lambda_step) = multiWavelengthLambdas;
	//TODO: provide Stokes calculation support
	if (lambda_start == lambda_end){
		//single wavelength calculation
		extJonesMain.reset(new ExtendedJones(materials, inAngles, lambda_start, lightSrcSpectrum, ifUseLambertian, false));
	}
	else{
		extJonesMain.reset(new ExtendedJones(materials, inAngles, lambda_start, lambda_end, lambda_step,
			lightSrcSpectrum, ifUseLambertian, false));
	}
}

std::vector<std::vector<std::pair<double, double> > > LCD1DMainBase::getIncidentAngles()const{
	std::vector<std::vector<std::pair<double, double> > > answer;
	std::vector<std::pair<double, double> > temp;
	std::pair<double, double> tp;
	for (auto& i : inAngles){
		temp.clear();
		for(auto& j : i){
			std::tie(tp.first, tp.second) = j;
			temp.push_back(tp);
		}
		answer.push_back(temp);
	}
	return answer;
}

LCD1DStaticMain::LCD1DStaticMain(double _lcLayerNum, double _dt, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing,
	double _voltStart, double _voltEnd, unsigned double _voltStep, double _maxIter, double _maxError)
	:LCD1DMainBase(_lcLayerNum, _dt, _lcParam, _rubbing), maxIter(_maxIter), maxError(_maxError){
		potentials.createStaticUpdatePolicy();
	}

LCD1DStaticMain::LCD1DStaticMain():LCD1DMainBase(){}

LCD::DOUBLEARRAY1D LCD1DStaticMain::getCalcVolts() const{
	return calcVolts;
}

std::vector<LCDOptics::TRANSRESULT> LCD1DStaticMain::getTransmissions()const{
	return transResults;
}

std::vector<DOUBLEARRAY2D> LCD1DStaticMain::getLCDirResults()const{
	return lcDirResult;
}

void LCD1DStaticMain::resetCalcVolts(double _voltStart, double _voltEnd, unsigned double _voltStep){
	calcVolts.clear();
	if (_voltEnd < _voltStart){
		double temp = _voltEnd;
		_voltEnd = _voltStart;
		_voltStart = temp;
	}

	for (double i = _voltStart; i <= _voltEnd + std::numeric_limits<double>::epsilon()*1000; i+=_voltStep){
		calcVolts.push_back(i);
	}
}

virtual void LCD1DStaticMain::calculate(){
	double residual = std::numeric_limits<double>::max();
	//empty the resut storage
	lcDirResult.clear();
	transResults.clear();
	for (auto volt: calcVolts){
		calculateOneVolt(volt);
		LCD1D::DIRVEC lcDir(lcDir->getDirectors());
		calc2X2OpticsOneSetLCDir(lcDir);
	}
}

double LCD1DStaticMain::calculateOneVolt(double _volt){
	double residual = std::numeric_limits<double>::max();
    //writeDirectors(lcDirector.getDirectors(), "_TNNoq0_"+toString(dc_volt)+"V_initialized_");
	unsigned long iternum = 0;
    while(residual >= maxError){
		if (iternum > maxIter){
			std::string msg;
			msg+="Static calculation does not converge when D.C. volt = ";
			msg+=toString(volt) + ".";
			throw std::runtime_error(msg.c_str());
		}
        potential.update(_volt);
        residual = lcDir->update();
		iternum++;
    };
	return residual;
}

void LCD1DStaticMain::calc2X2OpticsOneSetLCDir(LCD::DIRVEC lcDir){
	extJonesMain->resetLCDiretors(lcDir);
	extJonesMain->calculateExtendedJones();
	transResults.push_back(extJonesMain->getTransmissions());
}
