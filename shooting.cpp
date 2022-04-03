#include "shooting.hpp"

void ODE_solver::Second_order_step(double t)
{

    std::vector<double> k1 = f(yval,Parameters,t);
    
    for(size_t i=0;i<yval.size();i++) k1[i] *= step;
    for(size_t i=0; i< yval.size();i++) k1[i] = yval[i] + 0.5*k1[i];

    std::vector<double> k2 = f(k1,Parameters,t+0.5*step);

    for(size_t i=0; i< yval.size();i++) yval[i] += k2[i];

}

std::vector<double> ODE_solver::operator()()
{
    double t=start;
    while (t <= end)
    {
        Second_order_step(t);
        t+=step;
    }

    
    return yval;
}