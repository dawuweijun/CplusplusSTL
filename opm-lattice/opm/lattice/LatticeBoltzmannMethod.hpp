#ifndef LATTICEBOLTZMANNMETHOD_HEADER_INCLUDED
#define LATTICEBOLTZMANNMETHOD_HEADER_INCLUDED

#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/SimulatorState.hpp>
#include <vector>

class LatticeBoltzmannMethod {
public:
    LatticeBoltzmannMethod(const GridManager& grid, const LatticeBoltzmannModule& module, const FluidProperties& red, const FluidProperties& blue);
    ~LatticeBoltzmannMethod();
    void 
    step(const double dt, SimulatorState& state);
private:
    double externalForce_;
    double fluxForce_;
    double flux_;
    std::vector<double> gff_;
    std::vector<double> gfs_;
    const GridManager& grid_;
    const LatticeBoltzmannModule& module_;
    const FluidProperties& red_;
    const FluidProperties& blue_;
    void potential();
    void propagationBySwap(std::vector<std::vector<double>>& dist);
    void fcalcSc(std::vector<std::vector<double>>& f);
    double
    NipSc(const int flag, const double Rden, const double Bden, const double cxk, const double cyk, const double czk, const double wk, std::vector<double>& velocity);
    void collisionStepScBlue(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist);
    void collisionStepScRed(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist);
    void
    streamingSwap(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist);
};
#endif //LATTICE_HEADER_INCLUDED
