#include "Planet.h"
#include "PlanetMaker_Base.hpp"

extern "C" {
Planet Initialize(ModelType modtype)
{
    PlanetMaker_Base::Type Innertype;
    switch (modtype)
    {
    case ADIABATIC:
        Innertype = PlanetMaker_Base::Type::ADIABATIC;
        break;
    
    default:
        break;
    }
    PlanetMaker_Base* NewInst = new PlanetMaker_Base;
    *NewInst = PlanetMaker_Base::buildInstance(Innertype);
    return (PlanetMaker_Base*) NewInst;
}

void Delete(Planet planet)
{
    delete (PlanetMaker_Base*) planet;
}

DBtype BuildEnvelope(Planet planet, DBtype core, DBtype rho_neb, DBtype c_s)
{
    return ((PlanetMaker_Base*) planet)->calculateM_env(core,rho_neb,c_s);
}
 
}