#include <iostream>
#include "shooting.hpp"

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

    cout << iterator(1.5,1,2);

    return 0;
}