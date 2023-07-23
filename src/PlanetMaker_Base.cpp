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
    std::unique_ptr<PlanetBase> NewPlanet = nullptr;
    std::unique_ptr<InitPlanetBase> NewInitFunc = nullptr;
    std::unique_ptr<ScorePlanetBase> NewScore = nullptr;
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
    std::unique_ptr<Shooting_method> shoot{ new Shooting2(n,0,0,static_cast<Function2&>(*NewInitFunc),static_cast<Function&>(*NewPlanet),static_cast<MultiVariable&>(*NewScore))};

    std::unique_ptr<NewtonRaphson> NewSolver{ new NewtonRaphson{static_cast<AlgebraicFunction&>(*shoot)}};

    PlanetMaker_Base newOne(std::move(shoot),std::move(NewPlanet),std::move(NewInitFunc),std::move(NewScore),std::move(NewSolver));

    return newOne;
}



double PlanetMaker_Base::calculateM_env(Doub M_core, Doub rho_neb, Doub c_s)
{
    Doub M_env_guess{0};
    Doub R_core{std::pow((3*M_core / (4.* M_PI*5.5)),1./3.)};
    std::cerr<< "Start ready" << std::endl;
    InitPolitrop InitParams{rho_neb,0,0};
    PolitropParameters SetupParams;
    SetupParams.c_s=c_s;
    SetupParams.M_core=M_core;
    SetupParams.rho_neb=rho_neb*c_s*c_s;//actually P_neb... This is bad! 
    //std::shared_ptr<ParameterBase> InitPointer{&InitParams};
    Planet->setParams(&InitParams);
    ScoreFunc->M_core=M_core;
    

    //std::shared_ptr<ParameterBase> SetupPointer(&SetupParams);
    InitFunc->setup(&SetupParams);
    //Shooting_method shoot(2,GRAVI_CONST * (M_core+M_EARTH) / (c_s*c_s),R_core,*InitFunc,*Planet,*ScoreFunc);

    //NewtonRaphson NewSolver{shoot};
    

    CalcFunc->t2=R_core;
    CalcFunc->t1=GRAVI_CONST * (M_core+M_EARTH) / (c_s*c_s);
    std::cerr << "Init ready" << std::endl;

    //std::cerr << (*CalcFunc)(1*M_EARTH) << std::endl;
  /*  adiabatikus2 bolygo;
    bolygo.rho = rho_neb;
    adiabatikus_init tester_adiab{rho_neb*c_s*c_s,M_core,c_s};
    adiabatikus_score tester_adiab_score{M_core};

    Shooting2 tester_adiab_shoot{2,CalcFunc->t1,CalcFunc->t2,*InitFunc,*Planet,*ScoreFunc};
    
    NewtonRaphson tester_adiab_newton{tester_adiab_shoot};*/

    try
    {
        return (*Solver)(12 * M_EARTH, 1e-6 * M_EARTH , 15 * M_EARTH);

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

Envelope PlanetMaker_Base::BuildEnvelope(Doub M_core, Doub M_env, Doub rho_neb, Doub c_s)
{
    Doub R_core{std::pow((3*M_core / (4.* M_PI*5.5)),1./3.)};
    Doub R_out{GRAVI_CONST * (M_core+M_env) / (c_s*c_s)};
    std::vector<std::vector<double>> y_vec{std::vector<double>(1000),std::vector<double>(1000),std::vector<double>(1000)};
    PolitropParameters SetupParams;
    SetupParams.c_s=c_s;
    SetupParams.M_core=M_core;
    SetupParams.rho_neb=rho_neb*c_s*c_s;
    InitPolitrop InitParams{rho_neb,0,0};
    Planet->setParams(&InitParams);
    InitFunc->setup(&SetupParams);
    std::vector<double> y{rho_neb,M_env}, dy(2);
    (*InitFunc)(y,dy,&R_out);
    Second_order chosen_solver(static_cast<Function&>(*Planet),dy,R_out,R_core,(R_out-R_core)/1000);
    chosen_solver.solve_into_array(y_vec,ODE_solver::direction::BACKWARD);
    std::vector<double> radii = y_vec[0];
    std::vector<double> pressure = y_vec[1];
    std::vector<double> mass = y_vec[2];
    std::vector<double> rho(mass.size());
    std::vector<double> Temperature(mass.size());
    
    //rho = std::pow( P / K, gamma);
    std::transform(pressure.cbegin(),pressure.cend(),rho.begin(),
     [K=static_cast<adiabatikus2*>(Planet.get())->K,
     gam=static_cast<adiabatikus2*>(Planet.get())->gamma](double P){
        return std::pow(P/K,gam);});
    std::transform(pressure.cbegin(),pressure.cend(),Temperature.begin(),
     [K=static_cast<adiabatikus2*>(Planet.get())->K,
     gam=static_cast<adiabatikus2*>(Planet.get())->gamma](double P){
        return std::pow(P,1-gam)*std::pow(K,gam)*1.28/8.31432;});


    return Envelope{radii,mass,pressure,Temperature,rho};
}