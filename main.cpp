#include <iostream>
#include "shooting.hpp"

class Test_func : public Function{
    private:
    double a;
    std::vector<double> dy;

    public:
    std::vector<double> operator()(std::vector<double> y,std::vector<double> parameter, double t)
    {
        a = parameter[0];
        double b = parameter[1];
        dy[0] = a*y[0] + b;
        
    }

};


int main()
{


    return 0;
}