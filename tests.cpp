#include "shooting.hpp"
#include "Opacity.hpp"
#include <iostream>
#include <cmath>
#include<vector>
#include <iomanip>

class Test_func : public Function{
    private:
    double D;
    double m;

    public:
    Test_func(std::vector<double> params) : Function(), D(params[0]), m(params[1]){};

    void operator()(std::vector<double> y,std::vector<double> &dy, double t)
    {
        std::cerr << "Inner: " << y[0] << " " << y[1] << " "; // y.size() << " " << dy.size() << std::endl;
        dy[0] = y[1];
        dy[1] = - D/m*y[0];
        std::cerr << dy[0] << " " << dy[1] << std::endl;
    }

};

struct test_Algebraic : public AlgebraicFunction
{
    double operator()(double t)
    {
        return std::sin(t);
    }
};

struct Newton_test : public AlgebraicFunction
{
    double operator()(double x)
    {
        return x*x*x-2;
    }
};

int main()
{
    //Testing the fortran wrapper:

    int nT=101;
    char model[]="nrm";
    char top='h';
    char shape='s';
    bool ross=true;
    double rho0=1e-10;
    double rho1=1e-8;
    int nrho=4;
    double T0=10;
    double T1=9999;

    double* opactable= new double[nT*nrho];
    for(int i=0;i<nT*nrho;i++) opactable[i]=0;

    calculate_table_(model,&top,&shape,&ross,&T0,&T1,&nT,&rho0,&rho1,&nrho,opactable);

    std::ofstream test_opac{"opac_table.tab",std::ofstream::out | std::ofstream::trunc};

    test_opac << "T (K) rho ->" << std::scientific << std::setprecision(6);
    for(int i=0;i<nrho;i++) test_opac <<  "\t" << rho0*std::pow(std::exp(std::log(rho1/rho0)/(nrho-1)),static_cast<double>(i));
    test_opac << std::endl;

    for(int i=0;i<nT;i++)
    {
        test_opac << std::fixed << std::setprecision(4) << std::setw(12) <<
        T0*std::pow(std::exp(std::log(T1/T0)/(nT-1)),static_cast<double>(i))
          << std::setprecision(6);
        for(int j=0;j<nrho;j++)
        {
            test_opac << "\t" << opactable[i*nrho+j];
        }
        test_opac << std::endl;
    }
    test_opac.close();
    delete opactable;
}