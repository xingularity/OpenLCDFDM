#include <iostream>
#include "../../include/LCD_ExtendedJonesBase.hpp"

using namespace std;
using namespace LCDOptics;

int main(int argc, char const *argv[]) {
    LIGHTSPECTRUMDATA spec = LIGHTSPECTRUMDATA{{0.45, 0.7}, {0.55, 1.0}, {0.65, 0.83}};
    IAngles inAngles = std::vector<std::vector<Angle> >(1, std::vector<Angle>(2) );
    inAngles[0][0] = std::make_tuple(0.0,0.0);
    inAngles[0][1] = std::make_tuple(30.0/180.0*M_PI, 60.0/180.0*M_PI);
    NKData nkData;
    nkData[0.55] = COMPD(1.0, 0.0);
    NKoNKeData nkonkeData;
    nkonkeData[0.55] = std::make_tuple(COMPD(1.0,0.0), COMPD(2.0, 0.0));
    std::shared_ptr<Optical2X2OneLayer<ISOType> >isoLayers;
    isoLayers.reset(new Optical2X2OneLayer<ISOType>(1.0, nkData));
    std::shared_ptr<Optical2X2OneLayer<UniaxialType> >uniLayers;
    uniLayers.reset(new Optical2X2OneLayer<UniaxialType>(1.0, nkonkeData));
    std::shared_ptr<Optical2X2OneLayer<UniaxialType> >polarizer;
    polarizer.reset(new Optical2X2OneLayer<UniaxialType>(1.0, nkonkeData, OPT_POLARIZER));
    std:: cout << "print spectral efficiency" <<std::endl;
    for (auto i : SpectrumEfficiency){
        std::cout << i.first << ", " << i.second << std::endl;
    }
    MATERIALLAYERS2X2CONT layers;
    layers.push_back(isoLayers);
    layers.push_back(uniLayers);
    layers.push_back(polarizer);

    std::cout << "print material class:" << std::endl;
    for (auto i : layers)
        std::cout << i->opticalLayerKind() << std::endl;

    std::cout << "test single wavelength(0.55)" << std::endl;
    ExtendedJonesBase extJ1(layers, inAngles, 0.55, spec);
    std::cout << "LCLayer index:" << extJ1.getLCLayerIndex() << std::endl;
    std::cout << "print lambdas" << std::endl;
    for (auto i : extJ1.targetLambdas())
        std::cout << i << std::endl;
    std::cout << "print lightSourceSpectrum" << std::endl;
    for (auto i : extJ1.targetLightSrcSpectrum())
        std::cout << i << std::endl;
    std::cout << "print yBarOfLambda" << std::endl;
    for (auto i : extJ1.targetYBarOfLambda())
        std::cout << i << std::endl;
    std::cout << "Incident angles:" <<std::endl;
    for (auto i : extJ1.getIncidentAngles())
        for(auto j : i)
        std::cout << std::get<0>(j) << ", " << std::get<1>(j) << std::endl;

    std::cout << std::endl << "test multi-wavelength(0.3~0.85 step 0.46)" << std::endl;
    ExtendedJonesBase extJm(layers, inAngles, 0.3, 0.85, 0.046, spec);
    std::cout << "LCLayer index:" << extJm.getLCLayerIndex() << std::endl;
    std::cout << "print lambdas" << std::endl;
    for (auto i : extJm.targetLambdas())
        std::cout << i << std::endl;
    std::cout << "print lightSourceSpectrum" << std::endl;
    for (auto i : extJm.targetLightSrcSpectrum())
        std::cout << i << std::endl;
    std::cout << "print yBarOfLambda" << std::endl;
    for (auto i : extJm.targetYBarOfLambda())
        std::cout << i << std::endl;
    std::cout << "Incident angles:" <<std::endl;
    for (auto i : extJm.getIncidentAngles())
        for(auto j : i)
        std::cout << std::get<0>(j) << ", " << std::get<1>(j) << std::endl;

    return 0;
}
