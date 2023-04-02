#include <stdio.h>
#include "Planet.h"

#define M_EARTH  5.97216787e+27 //g

int main()
{
    Planet firstPlanet = Initialize(ADIABATIC);
    double rho_neb=1e-11;
    double M_core = 8* M_EARTH;
    double c_s = 148910;
    double M_env = BuildEnvelope(firstPlanet,M_core,rho_neb,c_s);
    printf("The calculated envelope of a %2f M_e planet is %4f M_e.\n ", M_core/M_EARTH, M_env/M_EARTH);
    Delete(firstPlanet);
    return 0;
}