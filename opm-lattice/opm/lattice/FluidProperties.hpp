#ifndef FLUIDPROPERTIES_HEADER_INCLUDED
#define FLUIDPROPERTIES_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include <opm/lattice/GridManager.hpp>

class FluidProperties{
public:
    FluidProperties(const GridManager& gird, const double rho, const double tau, const double mu, const double velmax);
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
};
#endif //FLUID_HEADER_INCLUDED
