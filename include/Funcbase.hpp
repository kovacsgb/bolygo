#ifndef FUNCBASE_DEFINED
#define FUNCBASE_DEFINED
#include <vector>
#include <memory>


typedef double Dvar;

struct ParameterBase
/** 
 * Placeholder in the Basic parameters struct. 
 * **/
{
    virtual ~ParameterBase() = default;
};
struct Function
{
     virtual void operator()(std::vector<double> y, std::vector<double> &dy, double t) = 0;
     virtual ~Function() {};
};

struct AlgebraicFunction
{
    virtual double operator()(double x) = 0;
    virtual ~AlgebraicFunction() {};
};

struct MultiVariable : public AlgebraicFunction
{
    virtual double operator()(double t, std::vector<double> x)= 0;

    virtual ~MultiVariable() {};
};

struct Function2 : public Function
{
    virtual void operator()(std::vector<double> y,std::vector<double> &dy, double t) {}
    virtual double operator()(std::vector<double> y,std::vector<double> &dy, double* t) {return 0;}
    virtual ~Function2() {};
};

struct PlanetBase : public Function
{

    virtual void setParams(ParameterBase *params)=0;
    virtual ParameterBase* getParams()=0;
    virtual ~PlanetBase() {};
};

struct InitPlanetBase : public Function2
{
    virtual void setup(ParameterBase* params)=0;
    virtual ~InitPlanetBase() {};
};

struct ScorePlanetBase : public MultiVariable
{
    Dvar M_core;

    ScorePlanetBase(Dvar x) : MultiVariable(), M_core(x) {}
    virtual ~ScorePlanetBase() = default;
};



#endif