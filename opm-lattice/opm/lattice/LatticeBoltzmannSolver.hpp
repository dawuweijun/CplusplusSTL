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
    ~LatticeBoltzmannSolver();
    void 
    step(const double dt, SimulatorState& state);
private:
    struct SolutionState {
        SolutionState(const int np);
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
    initVariables(const SimulatorState& state);
    void
    updateState(SimulatorState& state, const SolutionState& x);
    void 
    potential();
    void 
    propagationBySwap(SolutionState& x);
    void 
    fcalcSc(SolutionState& x);
    double
    NipSc(const int flag, const double Rden, const double Bden, const double cxk, const double cyk, const double czk, const double wk, std::vector<double>& velocity);
    void 
    collisionStepScBlue(SolutionState& x);
    void 
    collisionStepScRed(SolutionState& x);
    void
    streamingSwap(SolutionState& x);
    void
    massMomentumCalc(const SolutionState& x, SimulatorState& state);
    void
    pressureCalc(const SolutionState& x, SimulatorState& state);
};
#endif //LATTICEBOLTZMANNSOLVER_HEADER_INCLUDED
