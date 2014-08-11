#ifndef FLUIDPROPERTIES_HEADER_INCLUDED
#define FLUIDPROPERTIES_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <vector>
//class GridManager;
//class LatticeBoltzmannModule;
class FluidProperties{
public:
    FluidProperties(const GridManager& gird, const LatticeBoltzmannModule& module, const double, const double, const double, const double);
//    void init(const double x1, const double x2, const int loc, std::vector<double>& dist);
    const double rho() const { return  rho_; }
    const double tau() const { return tau_; }
    const double mu() const { return mu_; }
    const double velmax() const { return velmax_; }
    const double lambda() const { return -2.0 / (6*tau_*mu_ + 1); }
private:
    double rho_;
    double tau_;
    double mu_;
    double velmax_;
//    enum {Inside = 0, Outside = 1};
//    std::vector<double> velocity_;
//    std::vector<double> density_;
//    const GridManager& grid_;
//    const LatticeBoltzmannModule& module_;
//    std::vector<double> feq();
};
#endif //FLUID_HEADER_INCLUDED
