#include <stdio.h>

extern "C" {
extern void calculate_table_(char model[4],char* top, char* shape, bool* ross, double* T0,double* T1,int* NT,double* rho0,double* rho1,int* Nrho,double* opac);
/* This is the Cwrapper of the fortran calculate_table subroutine, called at initialization */
}
