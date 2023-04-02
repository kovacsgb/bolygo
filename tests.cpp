#include "shooting.hpp"
#include <iostream>
#include <cmath>
#include<vector>

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
