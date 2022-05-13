#include <iostream>
#include<iomanip>
#include<fstream>
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
    double const gamma = 1/(6./5.);
    double rho;


    void operator()(std::vector<double> y, std::vector<double> &dy, double t)
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
};

struct adiabatikus2 : public Function
{
    const double GRAVI_CONST = 6.674299999999999e-08;// in cgs -> cm3/gs2
    double const K = 3.85e12*0.5; // maybe and for the Sun, probably not the perfect.
    double gamma = 1/(4./3.);
    double rho;


    void operator()(std::vector<double> y, std::vector<double> &dy, double t)
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
};

struct adiabatikus_score : public MultiVariable
{
    double M_core;

    adiabatikus_score(double x2_) : M_core(x2_) {}

    double operator()(double t)
    {
        return 0;
    }

    double operator()(double t, std::vector<double> y)
    {
        return (double)((long double)M_core-(long double)y[1]);
    }

};

struct Function2 : public Function
{
    virtual void operator()(std::vector<double> y,std::vector<double> &dy, double t) {}
    virtual double operator()(std::vector<double> y,std::vector<double> &dy, double* t) {return 0;}
};

struct adiabatikus_init : public Function2
{
    const double GRAVI_CONST = 6.674299999999999e-08;
    double rho_neb;
    double M_core;
    double c_s;
    double R;

    adiabatikus_init(double rho_neb_, double M_core_,double c_s_) : rho_neb(rho_neb_),M_core(M_core_), c_s(c_s_) {}
    
    
    double operator()(std::vector<double> y,std::vector<double> &dy, double* t)
    {

        dy[0]=rho_neb; //actually P_neb
      //  std::cerr << "Szia uram!" << std::endl;
        dy[1]=M_core+y[1];
        R=GRAVI_CONST *(M_core+y[1])/(c_s*c_s);
        *t=R;
        return R;
    }
};




struct Shooting2 : public Shooting_method
{

    Function2& InitCalc2;
    std::ofstream log;

    Shooting2(int n, double t1_, double t2_, Function2& InitCalc_, Function& RHS_, MultiVariable& Score_) :
    Shooting_method(n,t1_,t2_,InitCalc_,RHS_,Score_), InitCalc2(InitCalc_), log("log"){}

    double operator()(double x)
{
    //--------------
    y[1] = x;
    dy[1] = 1000;
    //static_cast<double (&)(std::vector<double>,std::vector<double> &, double&)>(InitCalc.operator())(y,dy,t1);
    InitCalc2(y,dy,&t1);
    double h1 = (t2-t1)/1e3;
    log << t1 <<" " << t2 << " " << h1 << " " << dy[1] << " " << y[1] <<" " << x << std::endl;
    y=dy;
    Second_order integ(RHS,y,t1,t2,h1);
    try{
    y = integ();
    for (auto &&ye : y) std::cerr << ye;
    std::cerr << std::endl;
    for (auto &&ye : y) if (std::isnan(ye)) throw(0); 
    }
    catch(int e)
    {
        std::cerr << "One member of y is NAN! Write dump file.";
        std::ofstream dump{"dump"};
        integ.output.rdbuf(dump.rdbuf());
        integ(ODE_solver::direction::BACKWARD);
        throw("Not a Number error");
    }

    return Score(t2,y);
    //---------------
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
    double M_core = 10* M_EARTH;
    double r_core = std::pow((3*M_core / (4.* M_PI*5.5)),1./3.); //cm

    double M_env_guess = 1*M_EARTH;
    double M_tot = M_core+ M_env_guess;

    adiabatikus2 bolygo;
    bolygo.rho = rho_neb;
    double c_s = 148910;
    double R = bolygo.GRAVI_CONST * M_tot / (c_s*c_s);
    cerr << "#------------------" << endl;
    cerr << "M_tot=" << M_tot << " R=" << R << " r_core= " << r_core << " M_core= " << M_core << endl;

    
    //RK4 tryit(bolygo,{bolygo.K*std::pow(rho_neb,5./3.),R},M_tot,M_core,-(M_tot-M_core)/1e9);
    Second_order tryit_old(bolygo,{bolygo.K*std::pow(rho_neb,4.5/3.5),M_tot},R,r_core,-(R-r_core)/1e3);
    tryit_old(ODE_solver::direction::BACKWARD);

    adiabatikus_init tester_adiab{rho_neb*c_s*c_s,M_core,c_s};
    adiabatikus_score tester_adiab_score{M_core};

    Shooting_method tester_adiab_shoot{2,R,r_core,tester_adiab,bolygo,tester_adiab_score};
    
    NewtonRaphson tester_adiab_newton{tester_adiab_shoot};

    try
    {
        M_env_guess = tester_adiab_newton(1 * M_EARTH, 1e-6 * M_EARTH , 15 * M_EARTH);
    }
    catch(char const * e)
    {
        std::cerr << e << endl;
        return 0;
    }
    catch(string e)
    {
        std::cerr << e << endl;
        return 0;
    }
    M_tot=M_core+M_env_guess;
    R= bolygo.GRAVI_CONST * M_tot / (c_s*c_s);
    Second_order tryit(bolygo,{rho_neb*c_s*c_s,M_tot},R,r_core,-(R-r_core)/3e3);
    tryit(ODE_solver::direction::BACKWARD);


/*
    std::ofstream Fxtest{"Fxtest"};

    M_env_guess= -1 * M_EARTH + 1e-5 * M_EARTH;
    bolygo.gamma=1/(4./3.);
    while (M_env_guess < 2 * M_EARTH)
    {
        M_tot=M_core+M_env_guess;
        R= bolygo.GRAVI_CONST * M_tot / (c_s*c_s);
        //Second_order test_bolygo(bolygo,{bolygo.K*std::pow(rho_neb,1/bolygo.gamma),M_tot},R,r_core,-(R-r_core)/1e3);
        RK4 test_bolygo(bolygo,{rho_neb*c_s*c_s,M_tot},R,r_core,-(R-r_core)/1e3);
        auto y = test_bolygo(ODE_solver::direction::BACKWARD);
        Fxtest << M_env_guess << " " << tester_adiab_score(0, y) << endl;
        M_env_guess += 0.01 * M_EARTH;
    }
    
*/


    return 0;
}