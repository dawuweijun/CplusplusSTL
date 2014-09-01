#include <opm/lattice/FluidProperties.hpp>
#include <algorithm>
FluidProperties::FluidProperties(const LatticeBoltzmannModule& module, 
                                 const double rho, 
                                 const double tau, 
                                 const double mu, 
                                 const double velmax)
      :rho_(rho)
      ,tau_(tau)
      ,mu_(mu)
      ,velmax_(velmax)
      ,gff_(std::vector<double>(module.numDirection()))
      ,gfs_(std::vector<double>(module.numDirection()))
{
}

void
FluidProperties::setPotential(const double g0,
                              const double g1,
                              const double g2,
                              const double gs)
{
    double gff[19] = {  g1, g1, g1, g1, g1, g1,  
                        g2, g2, g2, g2, g2, g2,
                        g2, g2, g2, g2, g2, g2,
                        g0 };
    double gfs[19] = {  gs, gs, gs, gs, gs, gs,  
                        gs, gs, gs, gs, gs, gs,
                        gs, gs, gs, gs, gs, gs,  
                        g0 };
    
    std::copy(&gff[0], &gff[0]+19, gff_.begin());
    std::copy(&gfs[0], &gfs[0]+19, gfs_.begin());

}

