#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

LatticeBoltzmannModule::LatticeBoltzmannModule()
       : name_("D3Q19")
       , number_of_direction_(19)
       , weight_(std::vector<double>(19))
       , xVelocity_(std::vector<double>(19))
       , yVelocity_(std::vector<double>(19))
       , zVelocity_(std::vector<double>(19))
{
    const double w0 = 1.0 / 3.0;
    const double w1 = 1.0 / 36.0;
    const double w2 = 1.0 / 18.0;
    const double c0 = 0.0;
    const double c2 = 1.0;
    /* Default directional velocities */
    const double c_x[19] = {c2,-c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c2,-c2, c2,-c2, c0, c0, c0, c0, c0};
    const double c_y[19] = {c0, c0, c2,-c2, c0, c0, c2,-c2,-c2, c2, c0, c0, c0, c0, c2,-c2, c2,-c2, c0}; 
    const double c_z[19] = {c0, c0, c0, c0, c2,-c2, c0, c0, c0, c0, c2,-c2,-c2, c2, c2,-c2,-c2, c2, c0}; 
    /* Weight array */
    const double w[19] = {w2, w2, w2, w2, w2, w2, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w0};
    xVelocity_ = std::vector<double>(c_x, c_x + 19);
    yVelocity_ = std::vector<double>(c_y, c_y + 19);
    zVelocity_ = std::vector<double>(c_z, c_z + 19);
    weight_    = std::vector<double>(w  , w   + 19);
}

void
LatticeBoltzmannModule::initDistribution(const GridManager& grid, const FluidProperties& red, const FluidProperties& blue, SimulatorState& state) const
{
    const int ND = number_of_direction_;
    assert(ND == 19);
    std::vector<double> velocity(3, 0.0);
    
    std::cout << grid.NX() << std::endl;
    //set initial state.
    //create a tube.
    std::vector<int> innerIdx;
    std::vector<int> outerIdx;
    const int x1 = grid.NX() / 4;
    const int x2 = grid.NX()*3 / 4;
    const int xmax = std::max(x1, x2);
    const int xmin = std::min(x1, x2);
    std::cout << xmin << xmax<< std::endl;
    for (int z = 0; z < grid.NZ(); ++z) {
        for (int y = 0; y < grid.NY(); ++y) {
            for (int x = 0; x < grid.NX(); ++x) {
                if (x < xmin || x > xmax) {
                    std::cout << grid.index(x,y,z) <<"\n";
                    outerIdx.push_back(grid.index(x,y,z));
                } //else {
                if (x >= xmin && x <= xmax) {
                    innerIdx.push_back(grid.index(x,y,z));
                }
            }
        }
    }

    state.redDensity().resize(grid.dimension());
    state.blueDensity().resize(grid.dimension());
    std::vector<double> redDen(grid.dimension(), 0.0);
    std::vector<double> blueDen(grid.dimension(), 0.0);
    std::cout << "outer: " << outerIdx.size() << "inner: " << innerIdx.size() << "\n";
    for (auto i = 0; i < outerIdx.size(); ++i) {
        for (int k = 0; k < ND; ++k) {
            double uc = velocity[0]*xVelocity_[k] + velocity[1]*yVelocity_[k] + velocity[2]*zVelocity_[k];
    //        double uc = velocity[0] * cx[k] + velocity[1] * cy[k] + velocity[2] * cz[k];
            double u2 = std::pow(velocity[0], 2) + std::pow(velocity[1], 2) + std::pow(velocity[2],2);
            double feq = blue.rho()*weight_[k]*(1.0 + 3*uc + 4.5*std::pow(uc, 2) - 1.5*u2);
//            std::cout << outerIdx[i] <<"\n";
            blueDen[outerIdx[i]*ND + k] = feq; 
        }
    }
    for(auto i = 0; i < innerIdx.size(); ++i) {
        for (int k = 0; k < ND; ++k) {
            double uc = velocity[0]*xVelocity_[k] + velocity[1]*yVelocity_[k] + velocity[2]*zVelocity_[k];
            double u2 = std::pow(velocity[0], 2) + std::pow(velocity[1], 2) + std::pow(velocity[2],2);
            double feq = red.rho()*weight_[k]*(1.0 + 3*uc + 4.5*std::pow(uc, 2) - 1.5*u2);
            redDen[innerIdx[i]*ND + k] = feq;
        }
    }

    std::copy(redDen.begin(), redDen.end(), state.redDensity().begin());
    std::copy(blueDen.begin(), blueDen.end(), state.blueDensity().begin());
}
