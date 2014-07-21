#include "fluid.hpp"
#include <vector>


Fluid::Fluid(const Grid& grid, const Module& module, const double rho, const double tau, const double mu, const double velmax)
      :rho_(rho)
      ,tau_(tau)
      ,mu_(mu)
      ,velmax_(velmax)
      ,velocity_(std::vector<double>(grid.spaceDim()))
      ,density_(std::vector<double>(module.numDirection()))
      ,distribution_(std::vector<double>(grid.dimension(),std::vector<double>(module.numDirection())))
      ,gird_(grid)
      ,module_(module)
{
}

Fluid::~Fluid()
{
}

std::vector<double>
Fluid::feq()
{
    const int num = module_.numDirection();
    std::vector<double>& cx = module_.xVelocity();
    std::vector<double>& cy = module_.yVelocity();
    std::vector<double>& cz = module_.zVelocity();
    std::vector<double>& w = module_.weight();
    for (int k = 0; k < num; ++k) {
        double uc = velocity_[0] * cx[k] + velocity_[1] * cy[k] + velocity_[2] * cz[k];
        double u2 = std::pow(velocity_[0], 2) + std::pow(velocity_[1], 2) + std::pow(velocity_[2],2);
        density_[k] = rho_ * w[k] * (1.0 + 3 * uc + 4.5 * std::pow(uc, 2) - 1.5 * u2);
    }
    
    return density_;
}



Fluid::init(const double x1, const double x2, const Location loc)
{
    //set initial state.
    //create a tube.
    std::vector<int> idx;
    //create the index()
    for (int z = 0; z < grid.NZ(); ++z) {
        for (int y = 0; y < grid.NY(); ++y) {
            idx.push_back(grid.index(x1, y, z));
            idx.push_back(grid.index(x2, y, z));
        }
    }
    if (loc == Inside) {
        for(int i = 0; i < idx.size(); ++i) {
            distribution_[idx[i]] = feq();
        }
    }
    if (loc == Outside) {
        for (int i = 0; i < distribution_.size(); ++i) {
            for (auto elem : idx) {
                if (i != elem) {
                    distribution_[i] = feq();
                }
            }   
        }
    }
}
