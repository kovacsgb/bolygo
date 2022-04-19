#include <iostream>
#include<iomanip>
#include "shooting.hpp"


const double M_EARTH = 5.97216787e+27; // g

template<class T>
inline T toFour(T x)
{
    return x*x*x*x;
}


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


struct adiabatikus : public Function
{
    const double GRAVI_CONST = 6.674299999999999e-08;// in cgs -> cm3/gs2
    double const K = 3.85e15; // maybe and for the Sun, probably not the perfect.
    double const gamma = 1/(5./3.);
    double rho;


    void operator()(std::vector<double> y, std::vector<double> &dy, double t)
    {
        // dP/dM = - G M / 4 pi r^4
        // dr/dM = 1/4pi r^2 rho
        // P = K rho^gamma => rho = (P/K)^(1/gamma)
        double m = t;
        double P = y[0];
        double r = y[1];

        double dPdM = - ( GRAVI_CONST * m ) / (4 * M_PI * toFour(r));
        double drdM = 1 / (4* M_PI * r*r * rho);

        rho = std::pow( P / K, gamma);

        dy[0] = dPdM;
        dy[1] = drdM;
    }  
};



int main()
{
    using namespace std;
    vector<double> parameters = {1 ,1};
    vector<double> x0 = {1,0};
    //double t0, t1, ts;
    //t0=1; t1
    Test_func test(parameters);
    for (auto &&x : x0) cerr << x << " ";

    ODE_solver rugok(test,x0,0,10,0.01);

    rugok();
    cout << endl;
   
    Second_order midpoint(test,x0,0,10,0.01);
    midpoint();
    cout << endl;

    double t=30*M_PI;
    test_Algebraic tester;
    NumDeriv df(tester);
    while(t < 34* M_PI)
    {
        cout << t << " " << tester(t) << " " << df(t,tester(t)) << " " << df(t) - std::cos(t) << endl; 
    t+=0.01;
    }

    Newton_test tester2;
    NewtonRaphson iterator{tester2};

    cout << setprecision(10) << iterator(1.5,1,2);
    cout << endl << endl;

    test_init Load;
    test_score Score(2);
    Shooting_method shoot(2,0,2*M_PI,Load,test,Score);

    NewtonRaphson Shooter{shoot};
    cout << endl << "#--------------" << endl;

    double y0;
    try{

    y0 = Shooter(1,0.1,3);
    }
    catch(string ex)
    {
        cout << ex; 
    }

    cout << endl << "#--------" << endl;
    cout << y0;

    cout << endl << endl;

    //Here we try out bolygo

    double rho_neb=1e-11;
    double M_core = 1* M_EARTH;
    double r_core = std::pow((3*M_core / (4.* M_PI*5.5)),1./3.); //cm

    double M_env_guess = 800*M_EARTH;
    double M_tot = M_core+ M_env_guess;

    adiabatikus bolygo;
    bolygo.rho = rho_neb;
    double c_s = 148910;
    double R = bolygo.GRAVI_CONST * M_tot / (c_s*c_s);
    cerr << "#------------------" << endl;
    cerr << "M_tot=" << M_tot << " R=" << R << " r_core= " << r_core << " M_core= " << M_core << endl;

    Second_order tryit(bolygo,{bolygo.K*std::pow(rho_neb,5./3.),R},M_tot,M_core,-(M_tot-M_core)/1e4);
    tryit(ODE_solver::direction::BACKWARD);

    return 0;
}