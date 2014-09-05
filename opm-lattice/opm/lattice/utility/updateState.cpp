#include <opm/lattice/utility/updateState.hpp>
#include <vector>
#include <cassert>
#include <cmath>

void
updateState(SimulatorState& state, 
            const GridManager& grid, 
            const FluidProperties& red, 
            const FluidProperties& blue, 
            const LatticeBoltzmannModule& module)
{
    const int N = grid.dimension();
    const int ND = module.numDirection();
    const std::vector<int>& boundary = grid.boundary();
    for (int i = 0; i < N; ++i) {
        double tmp=0, tmp2 = 0;
        for (int k = 0; k < ND; ++k) {
            tmp += state.redDist()[i*ND + k];
            tmp2 += state.blueDist()[i*ND + k];
        }
        state.redDensity()[i] = tmp / (red.rho() + 1.0);
        state.blueDensity()[i] = tmp2 / (blue.rho() + 1.0);
        if (boundary[i] != 0) {
            state.redDensity()[i] = 0.0;
            state.blueDensity()[i] = 0.0;
        }
    }
}

