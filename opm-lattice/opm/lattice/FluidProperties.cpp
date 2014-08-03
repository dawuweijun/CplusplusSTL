#include <opm/lattice/FluidProperties.hpp>
#include <vector>


FluidProperties::FluidProperties(const GridManager& grid, const LatticeBoltzmannModule& module, const double rho, const double tau, const double mu, const double velmax)
      :rho_(rho)
      ,tau_(tau)
      ,mu_(mu)
      ,velmax_(velmax)
//      ,velocity_(std::vector<double>(grid.spaceDim()))
//      ,density_(std::vector<double>(module.numDirection()))
//      ,grid_(grid)
//      ,module_(module)
{
}


