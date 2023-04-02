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
    std::cout << "M_env =" << calculated_M_env/M_EARTH << std::endl;
    return 0;
}