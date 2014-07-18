#include "fluid.hpp"
#include <vector>

Fluid::Fluid(Grid& grid, double rho, double tau, double mu, double velmax)
{}

Fluid::Fluid(Grid& grid, Module& module, double rho, double tau, double mu, double velmax)
      :rho_(rho)
      ,tau_(tau)
      ,mu_(mu)
      ,velmax_(velmax)
      ,velocity_(grid.spaceDim())
      ,distribution_(module.numdirection(),std::vector<double>(grid.dimension()))
      ,gird_(grid)
      ,module_(module)
{
}

Fluid::~Fluid()
{}
std::vector<double>&
Fluid::feq()
{
    const int num = module_.numDirection();
    std::vector<double> pointDen(num);
    std::vector<double> cx = module_.xVelocity();
    std::vector<double> cy = module_.yVelocity();
    std::vector<double> cz = module_.zVelocity();
    std::vector<double> w = module_.weight();
    for (int k = 0; k < num; ++k) {
        double uc = velocity_[0] * cx[k] + velocity_[1] * cy[k] + velocity_[2] * cz[k];
        double u2 = std::pow(velocity_[0], 2) + std::pow(velocity_[1], 2) + std::pow(velocity_[2],2);
        pointDen[k] = rho_() * w[k] * (1.0 + 3 * uc + 4.5 * std::pow(uc, 2) - 1.5 * u2);
    }
    
    return pointDen;
}

Fluid::init()
{
    grid
}
