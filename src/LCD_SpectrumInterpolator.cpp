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
 *   * Neither the name of OpenLCDFDM nor
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

#include "LCD_OpticsBaseDef.hpp"

namespace LCDOptics{
    template <typename INTYPE>
    SpectrumInterpolator<INTYPE>::SpectrumInterpolator(const DOUBLEARRAY1D lambda_){static_assert(1, "This code shouldn't be instanitilized.");};

    template<>
    SpectrumInterpolator<LIGHTSPECTRUMDATA>::SpectrumInterpolator(){toleratedError = 1.0e-13;}

    template<>
    SpectrumInterpolator<LIGHTSPECTRUMDATA>::SpectrumInterpolator(const DOUBLEARRAY1D lambda_){lambdas = lambda_;toleratedError = 1.0e-13;}

    template<>
    SpectrumInterpolator<NKData>::SpectrumInterpolator(){toleratedError = 1.0e-13;};

    template<>
    SpectrumInterpolator<NKData>::SpectrumInterpolator(const DOUBLEARRAY1D lambda_){lambdas = lambda_;toleratedError = 1.0e-13;};

    template<>
    SpectrumInterpolator<NKoNKeData>::SpectrumInterpolator(){toleratedError = 1.0e-13;}

    template<>
    SpectrumInterpolator<NKoNKeData>::SpectrumInterpolator(const DOUBLEARRAY1D lambda_){lambdas = lambda_;toleratedError = 1.0e-13;};

    template<>
    SpectrumInterpolator<NK123Data>::SpectrumInterpolator(){toleratedError = 1.0e-13;};

    template<>
    SpectrumInterpolator<NK123Data>::SpectrumInterpolator(const DOUBLEARRAY1D lambda_){lambdas = lambda_;toleratedError = 1.0e-13;};

    template <typename T>
    std::tuple<double, double> SpectrumInterpolator<T>::searchLambdaRegion(const double lambda, const T& in){
        typename T::const_iterator iter = in.begin();
        typename T::const_iterator mapIterTail = in.end();
        --mapIterTail; //Go to the last element

        //check if the target wavelegnth is shorter than the shortest wavelength in nk distribution.
        //if true, returning a tuple with two values equal to shortest wavelength.
        if ((lambda - iter->first) <= 1.0*toleratedError) return std::make_tuple(iter->first, iter->first);

        //check if the target wavelegnth is longer than the longest wavelength in nk distribution.
        //if true, returning a tuple with two values equal to longest wavelength.
        if ((lambda - mapIterTail->first) >= -1.0*toleratedError) return std::make_tuple(mapIterTail->first, mapIterTail->first);

        //search otherwise.
        typename T::const_iterator tempRightIter = in.end();
        for (; iter != mapIterTail; iter++){
            //WATCH OUT! The left side iterator never reach --(in.end())
            tempRightIter = iter; tempRightIter++;
            if (((lambda - iter->first) > -1.0*toleratedError) && ((lambda - tempRightIter -> first) < 1.0*toleratedError)){
                return std::make_tuple(iter->first, tempRightIter->first);
            }
        }

        //should never reach here
        throw runtime_error("SpectrumInterpolator::searchLambdaRegion failed");
    }

    template <typename T>
    void SpectrumInterpolator<T>::interpolate(const T& nkData, T& wanted_nk){
        static_assert(1, "This code shouldn't be instanitilized.");
    }

    template<>
    void SpectrumInterpolator<LIGHTSPECTRUMDATA>::interpolate(const LIGHTSPECTRUMDATA& spectrum_in, LIGHTSPECTRUMDATA& targets){
        if(lambdas.size() == 0) throw runtime_error("lambdas has 0 elements while executing spectrumInterpolator<LIGHTSPECTRUMDATA>::interpolate().");
        if(spectrum_in.size() == 0) throw runtime_error("spectrum_in has 0 elements while executing SpectrumInterpolator<LIGHTSPECTRUMDATA>::interpolate().");

        LIGHTSPECTRUMDATA::const_iterator temp = spectrum_in.end();
        double minInputLambda = spectrum_in.begin()->first;
        double maxInputLambda = (--spectrum_in.end())->first;
        double ll,ul;//lower lambda, upper lambda

        targets.clear();
        for (DOUBLEARRAY1D::iterator iter = lambdas.begin(); iter != lambdas.end(); iter++){
            std::tie(ll,ul) = this->searchLambdaRegion(*iter, spectrum_in);
            //if the target wavelength <= (>=) shortest(longest) wavelength in spectrum_in
            if (ll == ul) targets[*iter] = spectrum_in.find(ll)->second;
            else targets[*iter] = (((spectrum_in.find(ul)->second))*(*iter - ll)+((spectrum_in.find(ll)->second))*(ul - *iter))/(ul-ll);
        }
    }

    template<>
    void SpectrumInterpolator<LIGHTSPECTRUMDATA>::interpolate(const LIGHTSPECTRUMDATA& spectrum_in, std::vector<LIGHTSPECTRUMDATA::value_type::second_type>& spectrum_out){
        LIGHTSPECTRUMDATA targets;
        this->interpolate(spectrum_in, targets);
        spectrum_out.clear();
        for (LIGHTSPECTRUMDATA::iterator iter = targets.begin(); iter != targets.end(); ++iter) spectrum_out.push_back(iter->second);
    }

    template<>
    LIGHTSPECTRUMDATA::value_type::second_type SpectrumInterpolator<LIGHTSPECTRUMDATA>::interpolate(const LIGHTSPECTRUMDATA& spectrum_in, double lambda){
        if(spectrum_in.size() == 0) throw runtime_error("spectrum_in has 0 elements while executing SpectrumInterpolator<LIGHTSPECTRUMDATA>::interpolate().");

        LIGHTSPECTRUMDATA::const_iterator temp = spectrum_in.end();
        double minInputLambda = spectrum_in.begin()->first;
        double maxInputLambda = (--spectrum_in.end())->first;
        double ll,ul;//lower lambda, upper lambda

        std::tie(ll,ul) = this->searchLambdaRegion(lambda, spectrum_in);
        //if the target wavelength <= (>=) shortest(longest) wavelength in spectrum_in
        if (ll == ul) return spectrum_in.find(ll)->second;
        else return (((spectrum_in.find(ul)->second))*(lambda - ll)+((spectrum_in.find(ll)->second))*(ul - lambda))/(ul-ll);
    }

    template<>
    NKData::value_type::second_type SpectrumInterpolator<NKData>::interpolate(const NKData& nkData, double lambda){
        if(nkData.size() == 0) throw runtime_error("nkData has 0 elements while executing SpectrumInterpolator<NKData>::interpolate().");

        NKData::const_iterator temp = nkData.end();
        double minInputLambda = nkData.begin()->first;
        double maxInputLambda = (--nkData.end())->first;
        double ll,ul;//lower lambda, upper lambda

        std::tie(ll,ul) = this->searchLambdaRegion(lambda, nkData);
        //if the target wavelength <= (>=) shortest(longest) wavelength in nkData
        if (ll == ul) return nkData.find(ll)->second;
        else return (((nkData.find(ul)->second))*(lambda - ll)+((nkData.find(ll)->second))*(ul - lambda))/(ul-ll);
    }

    template<>
    void SpectrumInterpolator<NKData>::interpolate(const NKData& nkData, NKData& targets){
        if(lambdas.size() == 0) throw runtime_error("lambdas has 0 elements while executing pectrumInterpolator<NKData>::interpolate().");
        if(nkData.size() == 0) throw runtime_error("nkData has 0 elements while executing SpectrumInterpolator<NKData>::interpolate().");

        NKData::const_iterator temp = nkData.end();
        double minInputLambda = nkData.begin()->first;
        double maxInputLambda = (--nkData.end())->first;
        double ll,ul;//lower lambda, upper lambda

        targets.clear();
        for (DOUBLEARRAY1D::iterator iter = lambdas.begin(); iter != lambdas.end(); iter++){
            std::tie(ll,ul) = this->searchLambdaRegion(*iter, nkData);
            //if the target wavelength <= (>=) shortest(longest) wavelength in nkData
            if (ll == ul) targets[*iter] = nkData.find(ll)->second;
            else targets[*iter] = (((nkData.find(ul)->second))*(*iter - ll)+((nkData.find(ll)->second))*(ul - *iter))/(ul-ll);
        }
    }

    template<>
    void SpectrumInterpolator<NKData>::interpolate(const NKData& nkData, std::vector<NKData::value_type::second_type>& wanted_nk){
        NKData targets;
        this->interpolate(nkData, targets);
        wanted_nk.clear();
        for (NKData::iterator iter = targets.begin(); iter != targets.end(); ++iter) wanted_nk.push_back(iter->second);
    }

    template<>
    NKoNKeData::value_type::second_type SpectrumInterpolator<NKoNKeData>::interpolate(const NKoNKeData& nkonkeData, double lambda){
        if(nkonkeData.size() == 0) throw runtime_error("nkonkeData has 0 elements while executing SpectrumInterpolator<NKoNKeData>::interpolate().");

        NKoNKeData::const_iterator temp = nkonkeData.end();
        double minInputLambda = nkonkeData.begin()->first;
        double maxInputLambda = (--nkonkeData.end())->first;
        double ll,ul;//lower lambda, upper lambda
        COMPD nko_l, nke_l;
        COMPD nko_u, nke_u;
        COMPD ans_nko, ans_nke;

        std::tie(ll,ul) = this->searchLambdaRegion(lambda, nkonkeData);
        //if the target wavelength <= (>=) shortest(longest) wavelength in nkonkeData
        if (ll == ul) return nkonkeData.find(ll)->second;
        else{
            std::tie(nko_l,nke_l) = nkonkeData.find(ll)->second;
            std::tie(nko_u,nke_u) = nkonkeData.find(ul)->second;
            ans_nko = ((nko_u)*(lambda - ll)+(nko_l)*(ul - lambda))/(ul-ll);
            ans_nke = ((nke_u)*(lambda - ll)+(nke_l)*(ul - lambda))/(ul-ll);
            return std::make_tuple(ans_nko, ans_nke);
        }
    }

    template<>
    void SpectrumInterpolator<NKoNKeData>::interpolate(const NKoNKeData& nkonkeData, NKoNKeData& targets){
        if(lambdas.size() == 0) throw runtime_error("lambdas has 0 elements while executing pectrumInterpolator<NKoNKeData>::interpolate().");
        if(nkonkeData.size() == 0) throw runtime_error("nkonkeData has 0 elements while executing SpectrumInterpolator<NKoNKeData>::interpolate().");

        NKoNKeData::const_iterator temp = nkonkeData.end();
        double minInputLambda = nkonkeData.begin()->first;
        double maxInputLambda = (--nkonkeData.end())->first;
        double ll,ul;//lower lambda, upper lambda
        COMPD nko_l, nke_l;
        COMPD nko_u, nke_u;
        COMPD ans_nko, ans_nke;

        targets.clear(); //a map container to put nk on target wavelengths
        for (DOUBLEARRAY1D::iterator iter = lambdas.begin(); iter != lambdas.end(); iter++){
            std::tie(ll,ul) = this->searchLambdaRegion(*iter, nkonkeData);
            //if the target wavelength <= (>=) shortest(longest) wavelength in nkonkeData
            if (ll == ul) targets[*iter] = nkonkeData.find(ll)->second;
            else{
                std::tie(nko_l,nke_l) = nkonkeData.find(ll)->second;
                std::tie(nko_u,nke_u) = nkonkeData.find(ul)->second;
                ans_nko = ((nko_u)*(*iter - ll)+(nko_l)*(ul - *iter))/(ul-ll);
                ans_nke = ((nke_u)*(*iter - ll)+(nke_l)*(ul - *iter))/(ul-ll);
                targets[*iter] = std::make_tuple(ans_nko, ans_nke);
            }
        }
    }

    template<>
    void SpectrumInterpolator<NKoNKeData>::interpolate(const NKoNKeData& nkonkeData, std::vector<NKoNKeData::value_type::second_type>& wanted_nk){
        NKoNKeData targets;
        this -> interpolate(nkonkeData, targets);
        wanted_nk.clear();
        for (NKoNKeData::iterator iter = targets.begin(); iter != targets.end(); ++iter) wanted_nk.push_back(iter->second);
    }

    template<>
    NK123Data::value_type::second_type SpectrumInterpolator<NK123Data>::interpolate(const NK123Data& nk123Data, double lambda){
        if(nk123Data.size() == 0) throw runtime_error("nk123Data has 0 elements while executing SpectrumInterpolator<NK123Data>::interpolate().");

        NK123Data::const_iterator temp = nk123Data.end();
        double minInputLambda = nk123Data.begin()->first;
        double maxInputLambda = (--nk123Data.end())->first;
        double ll,ul;//lower lambda, upper lambda
        COMPD nk1_l, nk2_l, nk3_l;
        COMPD nk1_u, nk2_u, nk3_u;
        COMPD ans_nk1, ans_nk2, ans_nk3;

        std::tie(ll,ul) = this->searchLambdaRegion(lambda, nk123Data);
        //if the target wavelength <= (>=) shortest(longest) wavelength in nk123Data
        if (ll == ul) return nk123Data.find(ll)->second;
        else{
            std::tie(nk1_l,nk2_l,nk3_l) = nk123Data.find(ll)->second;
            std::tie(nk1_u,nk2_u,nk3_u) = nk123Data.find(ul)->second;
            ans_nk1 = ((nk1_u)*(lambda - ll)+(nk1_l)*(ul - lambda))/(ul-ll);
            ans_nk2 = ((nk2_u)*(lambda - ll)+(nk2_l)*(ul - lambda))/(ul-ll);
            ans_nk3 = ((nk3_u)*(lambda - ll)+(nk3_l)*(ul - lambda))/(ul-ll);
            return std::make_tuple(ans_nk1, ans_nk2, ans_nk3);
        }
    }

    template<>
    void SpectrumInterpolator<NK123Data>::interpolate(const NK123Data& nk123Data, NK123Data& targets){
        if(lambdas.size() == 0) throw runtime_error("lambdas has 0 elements while executing spectrumInterpolator<NK123Data>::interpolate().");
        if(nk123Data.size() == 0) throw runtime_error("nk123Data has 0 elements while executing SpectrumInterpolator<NK123Data>::interpolate().");

        NK123Data::const_iterator temp = nk123Data.end();
        double minInputLambda = nk123Data.begin()->first;
        double maxInputLambda = (--nk123Data.end())->first;
        double ll,ul;//lower lambda, upper lambda
        COMPD nk1_l, nk2_l, nk3_l;
        COMPD nk1_u, nk2_u, nk3_u;
        COMPD ans_nk1, ans_nk2, ans_nk3;

        targets.clear(); //a map container to put nk on target wavelengths
        for (DOUBLEARRAY1D::iterator iter = lambdas.begin(); iter != lambdas.end(); iter++){
            std::tie(ll,ul) = this->searchLambdaRegion(*iter, nk123Data);
            //if the target wavelength <= (>=) shortest(longest) wavelength in nk123Data
            if (ll == ul) targets[*iter] = nk123Data.find(ll)->second;
            else{
                std::tie(nk1_l,nk2_l,nk3_l) = nk123Data.find(ll)->second;
                std::tie(nk1_u,nk2_u,nk3_u) = nk123Data.find(ul)->second;
                ans_nk1 = ((nk1_u)*(*iter - ll)+(nk1_l)*(ul - *iter))/(ul-ll);
                ans_nk2 = ((nk2_u)*(*iter - ll)+(nk2_l)*(ul - *iter))/(ul-ll);
                ans_nk3 = ((nk3_u)*(*iter - ll)+(nk3_l)*(ul - *iter))/(ul-ll);
                targets[*iter] = std::make_tuple(ans_nk1, ans_nk2, ans_nk3);
            }
        }
    }

    template<>
    void SpectrumInterpolator<NK123Data>::interpolate(const NK123Data& nk123Data, std::vector<NK123Data::value_type::second_type>& wanted_nk){
        NK123Data targets; //a map container to put nk on target wavelengths
        this-> interpolate(nk123Data, targets);
        wanted_nk.clear();
        for (NK123Data::iterator iter = targets.begin(); iter != targets.end(); ++iter) wanted_nk.push_back(iter->second);
    }
};
