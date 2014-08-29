#include <opm/lattice/FluidProperties.hpp>
FluidProperties::FluidProperties(const GridManager& grid, const double rho, const double tau, const double mu, const double velmax)
      :rho_(rho)
      ,tau_(tau)
      ,mu_(mu)
      ,velmax_(velmax)
{}


