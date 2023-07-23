#include "shooting.hpp"
#include <memory>
#include <algorithm>

typedef double Doub;


struct Parameters : public ParameterBase
{
    Doub rho_neb;
    Doub P_neb;
};

struct Envelope
{
    std::vector<Doub> Radius;
    std::vector<Doub> Mass;
    std::vector<Doub> Pressure;
    std::vector<Doub> Temperature;
    std::vector<Doub> Density;
};

class PlanetMaker_Base
{
    private:

    std::unique_ptr<Shooting_method> CalcFunc;
    std::unique_ptr<PlanetBase> Planet;
    std::unique_ptr<InitPlanetBase> InitFunc;
    std::unique_ptr<ScorePlanetBase> ScoreFunc;
    std::unique_ptr<NewtonRaphson> Solver;

    public:

    PlanetMaker_Base() : CalcFunc(nullptr), Planet(nullptr), InitFunc(nullptr), ScoreFunc(nullptr),Solver(nullptr) {}

    PlanetMaker_Base( std::unique_ptr<Shooting_method>&& method, std::unique_ptr<PlanetBase>&& planet, std::unique_ptr<InitPlanetBase>&& init, std::unique_ptr<ScorePlanetBase>&& score, std::unique_ptr<NewtonRaphson>&& solver):
       CalcFunc(std::move(method)), Planet(std::move(planet)), InitFunc(std::move(init)), ScoreFunc(std::move(score)), Solver(std::move(solver)) {}
    enum Type {ADIABATIC};

    static PlanetMaker_Base buildInstance(Type ModelType);

    double calculateM_env(Doub M_core, Doub rho_neb, Doub c_s);

    Envelope BuildEnvelope(Doub M_core, Doub M_env, Doub rho_neb, Doub c_s);
};