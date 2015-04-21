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

void LCD1DMainBase::setTFTPI(LCD1D::DielecParameters _tftpi={0.0, 0.0}){
	epsilonr->setTFTPI(_tftpi);
}

void LCD1DMainBase::setCFPI(LCD1D::DielecParameters _cfpi={0.0,0.0}){
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
