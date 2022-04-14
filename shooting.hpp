#ifndef SHOOTING_HEADER
#define SHOOTING_HEADER

#include <vector>
#include <iostream>



typedef double Dvar;
struct Function
{
     virtual void operator()(std::vector<double> y, std::vector<double> &dy, double t) = 0;
};


class ODE_solver{
/* Second order Runge-Kutta Solver
 * k1=h*f(x_n,y_n) 
 * k2 = h*f(x_n+1/2*h,y_n+1/2*k1)
 * y_n+1=y_n+k2
*/
    protected:
    std::vector<double> yval;
    std::vector<double> y0;
    std::vector<double> dy;
    //std::vector<double> Parameters; // Not necessary is part of Function
    double var;
    double start;
    double end;
    double step;
    Function& f;

    virtual void solv_step(double t);

    public:
    ODE_solver(); //contructor
    ODE_solver(Function& f_, std::vector<double> y0_,double start_,double end_,double step_) :
     yval(y0_.size()), y0(y0_), dy(y0_.size()), var(start_), start(start_), end(end_), step(step_), f(f_) {}

    std::vector<double> operator()();
    std::vector<double> operator()(std::vector<double> new_y0);

};

class Second_order : public ODE_solver
{
    std::vector<double> k1;
    void solv_step(double t);
    void Second_order_step(double t);

    public:
    Second_order(Function& f_, std::vector<double> y0_,double start_,double end_,double step_) : ODE_solver(f_,y0_,start_,end_,step_),  k1(y0_.size())
    {}

};

#endif