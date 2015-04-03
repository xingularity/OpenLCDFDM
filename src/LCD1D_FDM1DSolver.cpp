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

#include "LCD1D_FDM1DSolver.hpp"
#include "LCD_UsefulFuncs.hpp"
#include <stdexcept>

using namespace LCD;
using namespace LCD1D;

double Epsilon::calculateCapacitance() {
	if (epsr33.size() == 0) return 0.0;
	double answer=0.0;
	double piCap_inv=0.0;
	for (double epsr : epsr33){
		answer += dz/epsr;
	}
	lcCap = EPS0*1.0/answer;
	piCap_inv = tftpiThick/tftpi_epsr + cfpiThick/cfpi_epsr;
	if (std::abs(piCap_inv) < 1.0e-13)
		totalCap = lcCap;
	else
		totalCap = 1.0/(piCap_inv + 1.0/lcCap);
	return lcCap; //The unit is 1.0e-10F/cm^2
}

double Epsilon::getLCVoltRatio()const{
	double piCap{0.0};
	piCap = tftpiThick/tftpi_epsr + cfpiThick/cfpi_epsr;
	if (std::abs(piCap) < 1.0e-13)
		return 1.0;
	else{
		piCap = 1.0/piCap;
		return piCap / (piCap + lcCap);
	}
}

void Epsilon::updateEpsilonr(const double& epsr_para, const double& epsr_perp, const DIRVEC dirs){
	const double delta_epsr = epsr_para - epsr_perp;
	double theta1, theta2;
	if (dirs.extent(0) != (epsr33.size() + 1))
		throw std::runtime_error("Sizes of the director container and the eps33 container don't match.\n"
			+ std::string("dirs.extent(0) = ") + toString(dirs.extent(0)) + ", epsr33.size() = " + toString(epsr33.size()));
	for (int i = 0; i < epsr33.size(); ++i)
		epsr33[i] = epsr_perp + delta_epsr*std::pow(std::cos(std::acos(dirs(i)(2)) + std::acos(dirs(i+1)(2))), 2.0);
}

void LCDirector::resetConditions(const LCParamters _lcParam){
    lcParam = _lcParam;
    if (!lcUpdater)
        lcUpdater->changeLCParams(lcParam);
}

void LCDirector::resetConditions(const RubbingCondition _rubbing){
    rubbing = _rubbing;
    resetLCDirectors();
}

void LCDirector::createVectorFormUpdater(const Potential& _pot, double dt){
    lcUpdater.reset(new LCVecUpdater(_pot, *this, epsilonr, lcParam, cellgap/layerNum, dt));
}

LCSolverBase::LCSolverBase(const Potential& _pot, LCDirector& _lcDir, Epsilon& _epsilonr, LCParamters _lcParam
    , double _dz, double _dt):pot(_pot), epsilonr(_epsilonr), lcParam(_lcParam), lcDir(_lcDir), dz(_dz), dt(_dt){
        tempDir.resizeAndPreserve(lcDir.getSize());
}

LCVecUpdater::LCVecUpdater(const Potential& _pot, LCDirector& _lcDir, Epsilon& _epsilonr, LCParamters _lcParam, double _dz, double _dt):
    LCSolverBase(_pot, _lcDir, _epsilonr, _lcParam, _dz, _dt){
        temp.resizeAndPreserve(lcDir.getSize());
}

double LCVecUpdater::update(){
	//first, put the answer into temp and then swap.
	const DOUBLEARRAY1D& EFieldForLC = pot.getEfieldForLC();
	DIRVEC dirs = lcDir.lcDir; //use copy constructor
	double nx2=0.0, ny2=0.0, nz2=0.0, dnx=0.0, dny=0.0, dnz=0.0;

	//start to advancing directors
	for (int i=1; i < dirs.extent(0) - 1; i++){
		dnx=dirs(i+1)(0)-dirs(i-1)(0);
		dny=dirs(i+1)(1)-dirs(i-1)(1);
		dnz=dirs(i+1)(2)-dirs(i-1)(2);
		nx2=dirs(i)(0)*dirs(i)(0);
		ny2=dirs(i)(1)*dirs(i)(1);
		nz2=dirs(i)(2)*dirs(i)(2);

		//calculate variation of nx
		temp(i)(0)=((1.0-nz2)*k22+k33*nz2)*(dirs(i+1)(0)+dirs(i-1)(0)-
			2.0*dirs(i)(0))+q0*dz*k22*(dny)+0.5*dirs(i)(2)*(k33-k22)*(dnz*(dnx));
		temp(i)(0)/=(dz*dz*gamma);

		//calculate variation of ny
		temp(i)(1)=((1.0-nz2)*k22+k33*nz2)*(dirs(i+1)(1)+dirs(i-1)(1)-2.0*dirs(i)(1))
			-q0*dz*k22*(dnx)+0.5*dirs(i)(2)*(k33-k22)*(dnz*(dny));
		temp(i)(1)/=(dz*dz*gamma);
	
		//calculate variation of nz
		temp(i)(2)=(k11-nz2*k22+k33*nz2)*(dirs(i+1)(2)+dirs(i-1)(2)-2.0*dirs(i)(2))
			+EPS0*delta_epsr*EFieldForLC[i]*EFieldForLC[i]*dz*dz*dirs(i)(2) //divide dz later
			+0.25*dirs(i)(2)*(k33-k22)*(dnz*dnz-dnx*dnx-dny*dny);
		temp(i)(2)/=(dz*dz*gamma);
	}

	//update real directors
	for (int i=1;i < dirs.extent(0) - 1;i++){
		for (int j=0;j<3;j++)
			tempDir(i)(j) = dirs(i)(j) + dt*(temp(i)(j));
	}

	//normalize directors
	for (int i=1;i < dirs.extent(0) - 1;i++){
		double norm=0.0;
		for (int j=0;j<3;j++)
			norm+=tempDir(i)(j)*tempDir(i)(j);

		norm=std::sqrt(norm);
		for (int j=0;j<3;j++)
			tempDir(i)(j)/=norm;
	}

	//calculate the maximum residual
	double residual=0.0;
	double residual_temp = 0.0;
	for (int i=1;i < dirs.extent(0) - 1;i++){
		for (int j=0;j<3;j++){
			residual_temp=std::abs(tempDir(i)(j)-dirs(i)(j));
			residual = (residual_temp>residual)?residual_temp:residual;
		}
	}

	//put the advanced directors back to "dirs"
	blitz::swap(dirs, tempDir);
	epsilonr.updateEpsilonr(lcParam.epsr_para, lcParam.epsr_perp, dirs);
	return residual;
}

void Potential::update(double volt){potentialUpdater->update(volt);}

void Potential::createStaticUpdatePolicy(){
    potentialUpdater.reset(new PotentialSolversForStatic(potential, EFieldForLC, epsilonr, dz));
}
void Potential::createDynamicUpdatePolicy(std::shared_ptr<LCD::WaveformBase> _voltWavePtr){
    potentialUpdater.reset(new PotentialSolversForDynamic(potential, EFieldForLC, epsilonr, dz, _voltWavePtr));
}

PotentialCalculate::PotentialCalculate(DOUBLEARRAY1D& _pot, DOUBLEARRAY1D& _EFieldForLC, const Epsilon& _epsilons, double _dz):
	pot(_pot), EFieldForLC(_EFieldForLC), epsilonr(_epsilons), dz(_dz){

	matrixX3d.resize(_pot.size() - 2, 3);
	b.resize(_pot.size() - 2);
	for (int i = 0; i < b.size(); ++i)
		b(i) = 0.0;

	x.resize(_pot.size());
	for (auto& i : x) i = 0.0;
}

void PotentialCalculate::calculate(double volt){
	
	//apply B.C
	pot[0] = 0.0;
	pot.back() = volt;

	const DOUBLEARRAY1D& epsr33 = epsilonr.getLCEpsr33();
	
	int nGrid = pot.size();

	for (int i = 0; i < nGrid-2; i++)
		for (int j = 0; j < 3; j++)
			matrixX3d(i, j)=0.0;
	
	for (int i = 0; i < nGrid-2; i++){ b(i) = 0.0;}
	
	for (int i = 1; i < nGrid-3; i++){
		matrixX3d(i,0)=epsr33[i];
		matrixX3d(i,1)=(epsr33[i+1]+epsr33[i])*-1.0;
		matrixX3d(i,2)=epsr33[i+1];
	}

	//assign toppest and lowest row
	{
		matrixX3d(0,0)=(epsr33[1]+epsr33[0])*-1.0;
		matrixX3d(0,1)=epsr33[1];
		b(0)=-1.0*epsr33[0]*pot[0];//BC
	}

	{
		matrixX3d(nGrid-3, 0)=epsr33[nGrid-3];
		matrixX3d(nGrid-3, 1)=(epsr33[nGrid-2]+epsr33[nGrid-3])*-1.0;
		b(nGrid-3)=-1.0*epsr33[nGrid-2]*pot[nGrid-1];//BC
	}


	matrixX3d(1,1)-=matrixX3d(1,0)/matrixX3d(0,0)*matrixX3d(0,1);
	b(1)-=matrixX3d(1,0)/matrixX3d(0,0)*b(0);
	double temp=0.0;
	for (int i=1; i < nGrid-3;i++){
		temp=matrixX3d(i+1,0)/matrixX3d(i,1);
		matrixX3d(i+1,1)-=matrixX3d(i,2)*temp;
		b(i+1)-=b(i)*temp;
	}

	x[nGrid - 1] = volt;
	x[0] = 0.0;
	x[nGrid - 2] = b(nGrid-3)/matrixX3d(nGrid-3, 1);
	for (int i = nGrid-4; i>-1; i--)
		x[i+1]=(b(i)-matrixX3d(i,2)*x[i+2])/matrixX3d(i,1);

	//put the answers back.
	pot.swap(x);

	//update EFields
	for (size_t i = 1; i < pot.size() - 1; ++i)
		EFieldForLC[i] = (pot[i+1] - pot[i-1])/2.0/dz;
}

PotentialSolversForStatic::PotentialSolversForStatic(DOUBLEARRAY1D& _pot, DOUBLEARRAY1D& _EFieldForLC, const Epsilon& _epsilonr, double _dz)
	:PotentialCalculate(_pot, _EFieldForLC, _epsilonr, _dz){}

void PotentialSolversForStatic::update(double volt){
	//get potential B.C. on the LC layer
	volt*=epsilonr.getLCVoltRatio();
	this->calculate(volt);
}

PotentialSolversForDynamic::PotentialSolversForDynamic(DOUBLEARRAY1D& _pot, DOUBLEARRAY1D& _EFieldForLC, const Epsilon& _epsilons, double _dz, 
    std::shared_ptr<LCD::WaveformBase> _voltWavePtr)
	:PotentialCalculate(_pot, _EFieldForLC, _epsilons, _dz){
	voltWavePtr = _voltWavePtr;
	if (!voltWavePtr){
		std::cout << "No wave form applied, use 0.0 volt DC instead" << std::endl;
		voltWavePtr.reset(new DCWaveform(0.0));
	}
}

void PotentialSolversForDynamic::update(double t){
	this->calculate((*voltWavePtr)(t)*epsilonr.getLCVoltRatio());
}
