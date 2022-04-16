#include "shooting.hpp"


std::vector<double> ODE_solver::operator()(std::vector<double> new_y0)
{
    this->y0 = new_y0;
    return (*this)();
}




 std::vector<double> ODE_solver::operator()()
{
    double t=start;
    yval = y0;
    std::cerr << t << std::endl;
    while (t <= end)
    {
        this->solv_step(t);
        t+=step;
    }

    
    return yval;
}

void ODE_solver::solv_step(double t)
{
    f(yval,dy,t);
    std::cout << t << " ";
    for(auto&& y : yval)
    {
        std::cout << y << " ";
    }
    for(auto&& y: dy)
    {
        std::cout << y << " ";
    }
    std::cout << std::endl;

    for(size_t i=0;i<yval.size();i++) yval[i] += step*dy[i];
}


void Second_order::solv_step(double t)
{
    std::cout << t << " ";
    for(auto&& y : yval)
    {
        std::cout << y << " ";
    }
    for(auto&& y: dy)
    {
        std::cout << y << " ";
    }
    std::cout << std::endl;
    Second_order_step(t);
}

void Second_order::Second_order_step(double t)
{
    //Calculate k1
    f(yval,k1,t);
    //dy = k1;
    //for(size_t i=0;i<yval.size();i++) k1[i] *= step;
    for(size_t i=0; i< yval.size();i++) k1[i] = yval[i] + 0.5*step*k1[i];
    std::cerr << "y" << yval[0] << " " << yval[1] << " y+0.5*k1: ";
    for(auto &&k : k1) std::cerr << k << " ";

    //Calculate k2 := dy
    f(k1,dy,t+0.5*step);
    std::cerr << "k2= ";
    for(auto &&k : dy) std::cerr << k << " ";
    
    for(size_t i=0; i< yval.size();i++) yval[i] += step*dy[i];
    std::cerr << "new_y: ";
    for(auto &&k : yval) std::cerr << k << " ";
    std::cerr << t << " " << step << std::endl;

}

double NumDeriv::operator()(double x)
{
    double x_h = x; //use x+h
    volatile double temp = x_h; //need for roundoff error
    double h = temp != 0.0 ? prec * std::abs(temp) : prec ; //use x as curvature when calculate h
    x_h = temp + h;
    h = x_h - temp; //this is the machine precision trick

    double f_xph = f(x_h);
    double f_x = f(x);
    double derivative = (f_xph - f_x) / h;
    return derivative;
}

double NumDeriv::operator()(double x,double fx)
{
    double x_h = x; //use x+h
    volatile double temp = x_h; //need for roundoff error
    double h = temp != 0.0 ? prec * std::abs(temp) : prec ; //use x as curvature when calculate h
    x_h = temp + h;
    h = x_h - temp; //this is the machine precision trick

    double f_xph = f(x_h);
    double f_x = fx;
    double derivative = (f_xph - f_x) / h;
    return derivative;
}


double NewtonRaphson::operator()(double start, const double x1, const double x2)
{
    double xh,xl;
    double x = start;
    double f1=f(x1);
    double f2=f(x2);
    if ((f1 > 0.0 && f2 > 0.0) || (f1 <0.0 && f2 < 0.0)) throw ("Solution must be bracketed between x1,x2");
    if (f1 == 0.0) return x1;
    if (f2 == 0.0) return x2;
    if (f1 < 0.0)
    {
        xl=x1; xh = x2;
    }
    else
    {
        xh=x1; xl = x2;
    }
    double origrange = std::abs(x2-x1);
    double dx = origrange;
    double fx = f(start);
    double dfx = df(start,fx);

    for(int i=0;i<MAX_IT;i++)
    {
        //if out of range bisect:
        if ((((x-xh)*dfx-fx)*((x-xl)*dfx-fx)> 0) || (std::abs(2*fx) > std::abs(origrange*dfx)))
        {
            origrange = dx;
            dx = 0.5*(xh-xl);
            x = xl+dx;
            if(xl == x) return x;
        }
        else
        {
            origrange = dx;
            dx = fx/dfx;
            double temp = x;
            x -= dx;
            if(temp == x) return x;
        }
        if (std::abs(dx) < EPS) return x;
        fx =f(x);
        dfx=df(x,fx);
        if(fx < 0)
        {
            xl=x;
        }
        else
        {
            xh = x;
        }
    }

}