#ifndef PLANET
#define PLANET

#ifdef __cplusplus
extern "C"
{
#endif

typedef void* Planet;
typedef double DBtype;

typedef enum {ADIABATIC} ModelType;

Planet Initialize(ModelType modtype);

void Delete(Planet planet);

DBtype BuildEnvelope(Planet planet, DBtype core, DBtype rho_neb, DBtype c_s);


#ifdef __cplusplus
}
#endif


#endif