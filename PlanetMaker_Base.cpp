#include "PlanetMaker_Base.hpp"
#include <iostream>
#include <string>
#include<exception>
#include "shooting.hpp"
#include "politrop.hpp"


using namespace std;
const double M_EARTH = 5.97216787e+27;
const double GRAVI_CONST =6.674299999999999e-08;


class init_exception : public exception
{
    std::string funcname;
    int linenum;

    public:
    init_exception(const char* func, int line) : funcname(func), linenum(line) {}
    std::string what()
    {
        return std::string("Initialization error occured at ")+funcname+std::string(":")+std::to_string(linenum);
    }
};





PlanetMaker_Base PlanetMaker_Base::buildInstance(PlanetMaker_Base::Type ModelType)
{
    std::unique_ptr<Function> NewPlanet = nullptr;
    std::unique_ptr<Function2> NewInitFunc = nullptr;
    std::unique_ptr<MultiVariable> NewScore = nullptr;
    int n = 0;

    switch (ModelType)
    {
    case Type::ADIABATIC:
        NewPlanet   = std::unique_ptr<adiabatikus2>(new adiabatikus2);
        NewInitFunc = std::unique_ptr<adiabatikus_init>(new adiabatikus_init);
        NewScore    = std::unique_ptr<adiabatikus_score>(new adiabatikus_score);
        n=2;
        break;
    
    default:
        throw(init_exception(__func__,__LINE__));
        break;
    }
    std::unique_ptr<Shooting_method> shoot{ new Shooting_method(n,0,0,*NewInitFunc,*NewPlanet,*NewScore)};

    std::unique_ptr<NewtonRaphson> NewSolver{ new NewtonRaphson{*shoot}};

    PlanetMaker_Base newOne(std::move(shoot),std::move(NewPlanet),std::move(NewInitFunc),std::move(NewScore),std::move(NewSolver));

    return newOne;
}



double PlanetMaker_Base::calculateM_env(Doub M_core, Doub rho_neb, Doub c_s)
{
    Doub M_env_guess{0};
    Doub R_core{std::pow((3*M_core / (4.* M_PI*5.5)),1./3.)};
    std::cerr<< "Start ready" << std::endl;
    static_cast<adiabatikus2*>(Planet.get())->rho=rho_neb;
    static_cast<adiabatikus_score*>(ScoreFunc.get())->M_core=M_core;
    static_cast<adiabatikus_init*>(InitFunc.get())->setup(rho_neb*c_s*c_s,M_core,c_s);
    //Shooting_method shoot(2,GRAVI_CONST * (M_core+M_EARTH) / (c_s*c_s),R_core,*InitFunc,*Planet,*ScoreFunc);

    //NewtonRaphson NewSolver{shoot};
    

    CalcFunc->t2=R_core;
    CalcFunc->t1=GRAVI_CONST * (M_core+M_EARTH) / (c_s*c_s);
    std::cerr << "Init ready" << std::endl;
    //std::cerr << (*CalcFunc)(1*M_EARTH) << std::endl;
    adiabatikus2 bolygo;
    bolygo.rho = rho_neb;
    adiabatikus_init tester_adiab{rho_neb*c_s*c_s,M_core,c_s};
    adiabatikus_score tester_adiab_score{M_core};

    Shooting2 tester_adiab_shoot{2,CalcFunc->t1,CalcFunc->t2,tester_adiab,bolygo,tester_adiab_score};
    
    NewtonRaphson tester_adiab_newton{tester_adiab_shoot};

    try
    {
        return tester_adiab_newton(12 * M_EARTH, 1e-6 * M_EARTH , 15 * M_EARTH);
        //M_env_guess = (*Solver)(1 * M_EARTH, 1e-6 * M_EARTH , 15 * M_EARTH);
    }
    catch(char const * e)
    {
        std::cerr << e << endl;
        std::exit(0);
    }
    catch(string e)
    {
        std::cerr << e << endl;
        std::exit(0);
    }
    return M_env_guess;
}