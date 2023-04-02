#include "politrop.hpp"


void adiabatikus::operator()(std::vector<double> y, std::vector<double> &dy, double t)
    {
        // dP/dM = - G M / 4 pi r^450.5736	0^(1/gamma)
        double m = t;
        double P = y[0];
        double r = y[1];
        rho = std::pow( P / K, gamma);

        double dPdM = - ( GRAVI_CONST * m ) / (4 * M_PI * toFour(r));
        double drdM = 1 / (4* M_PI * r*r * rho);

        

        dy[0] = dPdM;
        dy[1] = drdM;
    }  


void adiabatikus2::operator()(std::vector<double> y, std::vector<double> &dy, double t)
    {
        // dP/dM = - G M / 4 pi r^450.5736	0^(1/gamma)
        double r = t;
        double P = y[0];
        double m = y[1];
        rho = std::pow( P / K, gamma);

        //double dPdM = - ( GRAVI_CONST * m ) / (4 * M_PI * toFour(r));
        //double drdM = 1 / (4* M_PI * r*r * rho);
        double dPdr = - GRAVI_CONST * m  * rho / (r*r);
        double dmdr = 4* M_PI * r*r* rho;

        dy[0] = dPdr;
        dy[1] = dmdr;
    }  

double adiabatikus_score::operator()(double t)
    {
        return 0;
    }

double adiabatikus_score::operator()(double t, std::vector<double> y)
    {
        return (double)((long double)M_core-(long double)y[1]);
    }


void adiabatikus_init::setup(double rho_neb_, double M_core_, double c_s_)
    {
        rho_neb=rho_neb_;
        M_core=M_core_;
        c_s=c_s_;
    }

double adiabatikus_init::operator()(std::vector<double> y,std::vector<double> &dy, double* t)
    {

        dy[0]=rho_neb; //actually P_neb
      //  std::cerr << "Szia uram!" << std::endl;
        dy[1]=M_core+y[1];
        R=GRAVI_CONST *(M_core+y[1])/(c_s*c_s);
        *t=R;
        return R;
    }
