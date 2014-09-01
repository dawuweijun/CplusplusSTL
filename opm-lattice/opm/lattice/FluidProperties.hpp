#ifndef FLUIDPROPERTIES_HEADER_INCLUDED
#define FLUIDPROPERTIES_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <vector>

class FluidProperties{
public:
    FluidProperties(const LatticeBoltzmannModule& module, 
                    const double rho, 
                    const double tau, 
                    const double mu, 
                    const double velmax);
    void setPotential(const double g0, 
                      const double g1, 
                      const double g2, 
                      const double gs);
    const double rho() const { return  rho_; }
    const double tau() const { return tau_; }
    const double mu() const { return mu_; }
    const double velmax() const { return velmax_; }
    const double lambda() const { return -2.0 / (6.0*tau_*mu_ + 1.0); }
    const std::vector<double> gff() const { return gff_; }
    const std::vector<double> gfs() const { return gfs_; }
private:
    double rho_;
    double tau_;
    double mu_;
    double velmax_;
    std::vector<double> gff_;
    std::vector<double> gfs_;    
};
#endif //FLUID_HEADER_INCLUDED
