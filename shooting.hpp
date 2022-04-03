#include <vector>

typedef double Dvar;
struct Function
{
     virtual std::vector<double> operator()(std::vector<double> y, std::vector<double> parameters, double t) = 0;
};


class ODE_solver{
/* Second order Runge-Kutta Solver
 * k1=h*f(x_n,y_n) 
 * k2 = h*f(x_n+1/2*h,y_n+1/2*k1)
 * y_n+1=y_n+k2
*/
    std::vector<double> yval;
    std::vector<double> y0;
    std::vector<double> Parameters;
    double var;
    double start;
    double end;
    double step;
    Function& f;

    ODE_solver(); //contructor

    std::vector<double> operator()();
    std::vector<double> operator()(std::vector<double> new_y0);

    void Second_order_step(double t);
};