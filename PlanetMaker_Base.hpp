#include "shooting.hpp"
#include <memory>

typedef double Doub;

struct Parameters
{
    Doub rho_neb;
    Doub P_neb;
};

class PlanetMaker_Base
{
    private:

    std::unique_ptr<Shooting_method> CalcFunc;
    std::unique_ptr<Function> Planet;
    std::unique_ptr<Function2> InitFunc;
    std::unique_ptr<MultiVariable> ScoreFunc;
    std::unique_ptr<NewtonRaphson> Solver;

    public:

    PlanetMaker_Base( std::unique_ptr<Shooting_method>&& method, std::unique_ptr<Function>&& planet, std::unique_ptr<Function2>&& init, std::unique_ptr<MultiVariable>&& score, std::unique_ptr<NewtonRaphson>&& solver):
       CalcFunc(std::move(method)), Planet(std::move(planet)), InitFunc(std::move(init)), ScoreFunc(std::move(score)), Solver(std::move(solver)) {}
    enum Type {ADIABATIC};

    static PlanetMaker_Base buildInstance(Type ModelType);

    double calculateM_env(Doub M_core, Doub rho_neb, Doub c_s);
};