#include <opm/lattice/utility/initState.hpp>
#include <vector>
#include <cassert>
#include <cmath>
void
initState(SimulatorState& state, 
          const GridManager& grid, 
          const FluidProperties& red, 
          const FluidProperties& blue, 
          const LatticeBoltzmannModule& module)
{
    const int ND = module.numDirection();
    assert(ND == 19);
    std::vector<double> velocity(3, 0.0);
    const std::vector<double>& cx = module.xVelocity();
    const std::vector<double>& cy = module.yVelocity();
    const std::vector<double>& cz = module.zVelocity();
    const std::vector<double>& w = module.weight();
    
    //set initial state.
    //create a tube.
    std::vector<int> innerIdx;
    std::vector<int> outerIdx;
    const int x1 = grid.NX() / 4;
    const int x2 = grid.NX()*3 / 4;
    const int xmax = std::max(x1, x2);
    const int xmin = std::min(x1, x2);
    for (int z = 0; z < grid.NZ(); ++z) {
        for (int y = 0; y < grid.NY(); ++y) {
            for (int x = 0; x < grid.NX(); ++x) {
                if (x < xmin || x > xmax) {
                    outerIdx.push_back(grid.index(x,y,z));
                } else {
                    innerIdx.push_back(grid.index(x,y,z));
                }
            }
        }
    }

    typedef std::vector<double>::size_type vec_size;
    for (vec_size i = 0; i < outerIdx.size(); ++i) {
        for (int k = 0; k < ND; ++k) {
            double uc = velocity[0] * cx[k] + velocity[1] * cy[k] + velocity[2] * cz[k];
            double u2 = std::pow(velocity[0], 2) + std::pow(velocity[1], 2) + std::pow(velocity[2],2);
            double feq = blue.rho()*w[k]*(1.0 + 3*uc + 4.5*std::pow(uc, 2) - 1.5*u2);
            state.blueDist()[outerIdx[i]*ND + k] = feq; 
        }
    }
    for(auto i = 0; i < innerIdx.size(); ++i) {
        for (int k = 0; k < ND; ++k) {
            double uc = velocity[0] * cx[k] + velocity[1] * cy[k] + velocity[2] * cz[k];
            double u2 = std::pow(velocity[0], 2) + std::pow(velocity[1], 2) + std::pow(velocity[2],2);
            double feq = red.rho()*w[k]*(1.0 + 3*uc + 4.5*std::pow(uc, 2) - 1.5*u2);
            state.redDist()[innerIdx[i]*ND + k] = feq;
        }
    }
}
