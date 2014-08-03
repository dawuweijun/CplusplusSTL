#ifndef LATTICEBOLTZMANNSOLVER_HEADER_INCLUDED
#define LATTICEBOLTZMANNSOLVER_HEADER_INCLUDED

#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/SimulatorState.hpp>
#include <vector>

class LatticeBoltzmannSolver {
public:
    LatticeBoltzmannSolver(const GridManager& grid, const LatticeBoltzmannModule& module, const FluidProperties& red, const FluidProperties& blue);
    void 
    step(const double dt, SimulatorState& x);
private:
    struct SolutionState {
        SolutionState(const int size);
        std::vector<double> red;
        std::vector<double> blue;
    };
    double externalForce_;
    double fluxForce_;
    double flux_;
    std::vector<double> gff_;
    std::vector<double> gfs_;
    const GridManager& grid_;
    const LatticeBoltzmannModule& module_;
    const FluidProperties& red_;
    const FluidProperties& blue_;
    SolutionState
    initVariables();
    void
    initDistribution(const double rho, const double x1, const double x2, const int loc, std::vector<double>& dist);
    void
    updateState(const SolutionState& state, SimulatorState& x);
    void 
    potential();
    void 
    propagationBySwap(std::vector<double>& dist);
    void 
    fcalcSc(std::vector<double>& f);
    double
    NipSc(const int flag, const double Rden, const double Bden, const double cxk, const double cyk, const double czk, const double wk, std::vector<double>& velocity);
    void 
    collisionStepScBlue(SolutionState& state);
    void 
    collisionStepScRed(SolutionState& state);
    void
    streamingSwap(SolutionState& state);
    void
    massMomentumCalc(const SolutionState& state, SimulatorState& x);
    void
    pressureCalc(const SolutionState& state, SimulatorState& x);
};
#endif //LATTICEBOLTZMANNSOLVER_HEADER_INCLUDED
