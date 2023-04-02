#include "Funcbase.hpp"
#include <cmath>

template<class T>
inline T toFour(T x)
{
    return x*x*x*x;
}


struct adiabatikus : public Function
{
    const double GRAVI_CONST = 6.674299999999999e-08;// in cgs -> cm3/gs2
    double const K = 3.85e15; // maybe and for the Sun, probably not the perfect.
    double const gamma = 1/(6./5.);
    double rho;


    void operator()(std::vector<double> y, std::vector<double> &dy, double t);
 };

struct adiabatikus2 : public Function
{
    const double GRAVI_CONST = 6.674299999999999e-08;// in cgs -> cm3/gs2
    double const K = 3.85e12*0.5; // maybe and for the Sun, probably not the perfect.
    double gamma = 1/(4./3.);
    double rho;


    void operator()(std::vector<double> y, std::vector<double> &dy, double t);
};

struct adiabatikus_score : public MultiVariable
{
    double M_core;

    adiabatikus_score(double x2_) : M_core(x2_) {}
    adiabatikus_score() : M_core(0) {}

    double operator()(double t);
    double operator()(double t, std::vector<double> y);

};



struct adiabatikus_init : public Function2
{
    const double GRAVI_CONST = 6.674299999999999e-08;
    double rho_neb;
    double M_core;
    double c_s;
    double R;

    adiabatikus_init(double rho_neb_, double M_core_,double c_s_) : rho_neb(rho_neb_),M_core(M_core_), c_s(c_s_) {}
    adiabatikus_init() {}
    
    void setup(double rho_neb_, double M_core_, double c_s_);

    double operator()(std::vector<double> y,std::vector<double> &dy, double* t);
};

