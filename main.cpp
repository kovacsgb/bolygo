#include <iostream>
#include<iomanip>
#include<fstream>
#include "shooting.hpp"
#include "PlanetMaker_Base.hpp"


const double M_EARTH = 5.97216787e+27; // g

int main()
{
  
    PlanetMaker_Base instance=PlanetMaker_Base::buildInstance(PlanetMaker_Base::Type::ADIABATIC);
    double rho_neb=1e-11;
    double M_core = 8* M_EARTH;
    double c_s = 148910;
    double calculated_M_env = instance.calculateM_env(M_core,rho_neb,c_s);
    Envelope full_envelope  = instance.BuildEnvelope(M_core,calculated_M_env,rho_neb,c_s);

    std::ofstream envOut{"Envelope.dat",std::ofstream::out | std::ofstream::trunc};
    if(envOut.is_open()) std::cout << "I am open!" << std::endl;
    envOut << std::scientific << std::setprecision(6) << std::setw(10);
    for(size_t i=0;i<full_envelope.Radius.size();i++)
    {
        envOut << full_envelope.Radius[i] << " "
               << full_envelope.Mass[i]   << " "
               << full_envelope.Pressure[i] << " "
               << full_envelope.Temperature[i] << " "
               << full_envelope.Density[i] << " "
               << std::endl;
    }
    envOut.close();
    std::cout << "M_env =" << calculated_M_env/M_EARTH << std::endl;
    return 0;
}