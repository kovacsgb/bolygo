#ifndef SHOOTING_HEADER
#define SHOOTING_HEADER

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <limits>
#include <functional>
#include "Funcbase.hpp"



class ODE_solver{
/* Second order Runge-Kutta Solver
 * k1=h*f(x_n,y_n) 
 * k2 = h*f(x_n+1/2*h,y_n+1/2*k1)
 * y_n+1=y_n+k2
*/
    private:
    std::vector<double> Forward;

    protected:
    std::vector<double> yval;
    std::vector<double> y0;
    std::vector<double> dy;
    double var;
    double start;
    double end;
    double step;
    Function& f;
    std::function<std::vector<double>()> solving;

    virtual void solv_step(double t);

    public:
    enum direction { FORWARD, BACKWARD};
    std::ostream& output;
    ODE_solver(); //contructor
    ODE_solver(Function& f_, std::vector<double> y0_,double start_,double end_,double step_);

    std::vector<double> operator()();
    std::vector<double> operator()(std::vector<double> new_y0);
    std::vector<double> operator()(direction where);
    void solve_into_array(std::vector<std::vector<double>> &results, direction where);
    

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

class RK4 : public ODE_solver
{
    std::vector<double> k1;
    std::vector<double> k2;
    std::vector<double> k3;
    std::vector<double> k4;
    std::vector<double> temp;

    void solv_step(double t);
    void RK4_step(double t);

    public:
    RK4(Function& f_, std::vector<double> y0_,double start_,double end_,double step_) : 
    ODE_solver(f_,y0_,start_,end_,step_),  k1(y0_.size()), k2(y0_.size()), k3(y0_.size()),
    k4(y0_.size()), temp(y0_.size())
    {}


};

struct NumDeriv
/*
 * Functor calculating numerical derivative.
 * Error of the derivative calculation is sqrt(epsilon_m)*x_c ~ 1e-8 * x
 * Better error can be achieved by more function evaulation;
 */
{
    AlgebraicFunction& f;
    const double prec;

    NumDeriv(AlgebraicFunction& f_) : f(f_), prec(calc_prec<double>()) {}

    double operator()(double x);
    //Retun df/dx at x
    double operator()(double x, double fx);
    //Return df/dx with single function evaulation.

    private:
    template<class T>
    double calc_prec()
    {
        double machine_precision = std::numeric_limits<T>::epsilon();
        return std::sqrt(machine_precision);
    }

};



class NewtonRaphson
/*
    Newton-Raphson iterator with bisection method as safeguard.
    In the original problem we have only 1 variable M_env which is not known in the outer edge of the planet, since R is dependent on M_tot.
    In this case we only need 1D newton iterator.
*/
{
    private:
    const int MAX_IT = 200;
    const double EPS = std::sqrt(std::numeric_limits<double>::epsilon()); //maybe this will be better //I am not sure this will be good, if error in df is sqrt(EPS).
    AlgebraicFunction& f;
    NumDeriv df;



    public:
    NewtonRaphson(AlgebraicFunction& f_) : f(f_), df(f_) {}

    double operator()(double start,const double x1, const double x2);

};

struct Shooting_method : public AlgebraicFunction
/* Original base class for the Shooting method, uses forward calculation. */
{
    double t1,t2;
    std::vector<double> init_vals;
    std::vector<double> y;
    std::vector<double> dy;
    Function& RHS;
    Function& InitCalc;
    MultiVariable& Score;

    Shooting_method(int n, double t1_, double t2_, Function& InitCalc_, Function& RHS_, MultiVariable& Score_) :
     t1(t1_),t2(t2_), init_vals(n), y(n), dy(n), RHS(RHS_), InitCalc(InitCalc_), Score(Score_) {}

    virtual double operator()(double x);
};

struct Shooting2 : public Shooting_method
/* Backward calculating shooting method. This works well with polytrope envelope.
*/
{

    Function2& InitCalc2;
    std::ofstream log;

    Shooting2(int n, double t1_, double t2_, Function2& InitCalc_, Function& RHS_, MultiVariable& Score_) :
    Shooting_method(n,t1_,t2_,InitCalc_,RHS_,Score_), InitCalc2(InitCalc_), log("log"){}

    double operator()(double x) override;


};

struct test_score : public MultiVariable
{
    double x2;

    test_score(double x2_) : x2(x2_) {}

    double operator()(double t)
    {
        return 0;
    }

    double operator()(double t, std::vector<double> y)
    {
        return x2-y[0];
    }

};

struct test_init : public Function
{
    void operator()(std::vector<double> y,std::vector<double> &dy, double t)
    {
        dy[0]=y[0];
        dy[1]=0;
    }
};


#endif