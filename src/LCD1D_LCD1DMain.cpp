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

LCD1DMainBase::LCD1DMainBase(double _lcLayerNum, LCD1D::LCParamters _lcParam, LCD1D::RubbingCondition _rubbing){
	//TODO
}

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

void LCD1DMainBase::addOpticalGlassLayer(double _thick, std::vector<double> _spectrumLambdas, std::vector<std::complex<double> > _nkSpectrum){
	addOpticalIsotropicLayer(_thick, _spectrumLambdas, _nkSpectrum, LCDOptics::OPT_GLASS);
}

void LCD1DMainBase::addOpticalIsotropicLayer(double _thick, std::vector<double> _spectrumLambdas,
	std::vector<std::complex<double> > _nkSpectrum, LCDOptics::OpticalMaterialClass _class){

	if (_spectrumLambdas.size() != _nkSpectrum.size())
		throw runtime_error("Size of the wavelengths don't corresspond to the size of nk values");

	LCDOptics::NKData nkSpectrum;
	for (int i = 0; i < _spectrumLambdas.size(); ++i)
		nkSpectrum[_spectrumLambdas[i]] = _nkSpectrum[i];

	materials.push_back(LCDOptics::Optical2x2IsoPtr(
		new LCDOptics::Optical2x2IsoPtr::element_type(thick, nkSpectrum, _class)));
}

void LCD1DMainBase::addOpticalUnaixialLayer(double _thick, std::vector<double> _spectrumLambdas,
	std::vector<std::complex<double> > _nokoSpectrum, std::vector<std::complex<double> > _nekeSpectrum,
	LCDOptics::OpticalMaterialClass _class){

	std::vector<std::complex<double> > _nkSpectrum, LCDOptics::OpticalMaterialClass _class){
	if (_spectrumLambdas.size() != _nokoSpectrum.size())
		throw runtime_error("Size of the wavelengths don't corresspond to the size of nk values");
	else if (_spectrumLambdas.size() != _nekeSpectrum.size())
		throw runtime_error("Size of the wavelengths don't corresspond to the size of nk values");

	LCDOptics::NKoNKeData nkSpectrum;
	for (int i = 0; i < _spectrumLambdas.size(); ++i){
		COMPD nko(_nokoSpectrum[i].real(), -1.0*_nokoSpectrum[i].imag());
		COMPD nke(_nekeSpectrum[i].real(), -1.0*_nekeSpectrum[i].imag());
		//make sure the k value us negative
		nko.imag() =
		nkSpectrum[_spectrumLambdas[i]] = std::make_tuple(nko, nke);
	}
	
	materials.push_back(LCDOptics::Optical2x2IsoPtr(
		new LCDOptics::Optical2x2IsoPtr::element_type(thick, nkSpectrum, _class)));
}

void LCD1DMainBase::addOpticalPolarizer(double _thick, std::vector<double> _spectrumLambdas, std::vector<std::complex<double> > _nokoSpectrum, std::vector<std::complex<double> > _nekeSpectrum){
    addOpticalUnaixialLayer(_thick, _spectrumLambdas, _nokoSpectrum, _nekeSpectrum, LCDOptics::OPT_POLARIZER);
}

void LCD1DMainBase::addOpticalLC(double _thick, std::vector<double> _spectrumLambdas, std::vector<std::complex<double> > _nokoSpectrum, std::vector<std::complex<double> > _nekeSpectrum){
    addOpticalUnaixialLayer(_thick, _spectrumLambdas, _nokoSpectrum, _nekeSpectrum, LCDOptics::OPT_LCMATERIAL);
}

void LCD1DMainBase::setOpticalIncidentAngles(unsigned double _intervalDegree){
    std::vector<double> thetas(1,0.0);
    std::vector<double> phis(1, 0.0);
    
    while(thetas.back()+_intervalDegree <= (360+1.0e-13)){
        thetas.push_back(thetas.back() + _intervalDegree);
    };
    
    while(phis.back()+_intervalDegree <= (80+1.0e-13)){
        phis.push_back(phis.back() + _intervalDegree);
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

void LCD1DMainBase::setOpticalWavelength(double _lambda_start, double _lambda_end, double _lambda_step){
    multiWavelengthLambdas = std::make_tuple(_lambda_start, _lambda_end, _lambda_step);
}
void LCD1DMainBase::setOpticalWavelength(double _lambda){
    multiWavelengthLambdas = std::make_tuple(_lambda, _lambda, 0);
}

void LCD1DMainBase::disableOpticalCalculation(){
    multiWavelengthLambdas = std::make_tuple(0.0,0.0,0.0);
}

void LCD1DMainBase::setOpticalSourceSpectrum(std::vector<double> _lambdas, std::vector<double> _powers){
    if (_lambdas.size() != _power.size())
        throw std::runtime_error("_lambdas.size() != _power.size() in LCD1DMainBase::setOpticalSourceSpectrum");
    lightSrcSpectrum.clear();
    for (int i = 0; i < _lambdas.size(); ++i)
        lightSrcSpectrum[_lambdas[i]] = _powers[i];
}

void LCD1DMainBase::calculate(){
    
}

void LCD1DMainBase::createExtendedJones(){

}
