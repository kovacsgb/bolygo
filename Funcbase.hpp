#ifndef FUNCBASE_DEFINED
#define FUNCBASE_DEFINED
#include <vector>


typedef double Dvar;

struct Function
{
     virtual void operator()(std::vector<double> y, std::vector<double> &dy, double t) = 0;
};

struct AlgebraicFunction
{
    virtual double operator()(double x) = 0;
};

struct MultiVariable : public AlgebraicFunction
{
    virtual double operator()(double t, std::vector<double> x)= 0;
};

struct Function2 : public Function
{
    virtual void operator()(std::vector<double> y,std::vector<double> &dy, double t) {}
    virtual double operator()(std::vector<double> y,std::vector<double> &dy, double* t) {return 0;}
};

#endif