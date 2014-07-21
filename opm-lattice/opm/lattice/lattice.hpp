#ifndef LATTICE_HEADER_INCLUDED
#define LATTICE_HEADER_INCLUDED

#include <opm/lattice/grid.hpp>
#include <opm/lattice/module.hpp>
#include <opm/lattice/fluid.hpp>
#include <vector>

class Lattice {
public:
    Lattice(const Grid& grid, const Module& module, const Fluid& red, const Fluid& blue);
    ~Lattice();
    void 
    step(const Grid& grid, Fluid& red, Fluid& blue);
private:
    double externalForce_;
    double fluxForce_;
    double flux_;
    std::vector<double> gff_;
    std::vector<double> gfs_;
    const Grid& grid_;
    const Module& module_;
    const Fluid& red_;
    const Fluid& blue_;
    void potential();
    void propagationBySwap(std::vector<std::vector<double>>& dist);
    void fcalcSc(std::vector<std::vector<double>>& f);
    double
    NipSc(const int flag, const double Rden, const double Bden, const double cxk, const double cyk, const double czk, const double wk, std::vector<double>& velocity);
    void collisionStepScBlue(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist);
    void collisionStepScRed(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist);
    void
    streamingSwap(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist);
    void massMomentumCalc(const std::vector<Fluid>& fluid, const Grid& grid);
    void pressureCalc(const std::vector<Fluid>& fluid, const Grid& grid);
};
#endif //LATTICE_HEADER_INCLUDED
