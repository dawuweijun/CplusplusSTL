#include <opm/lattice/SimulatorState.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>

void SimulatorState::init(const GridManager& grid, const LatticeBoltzmannModule& module)
{
    num_phases_ = 2;
    press_.resize(grid.dimension(), 0.0);
    vel_.resize(num_phases_ * grid.dimension(), 0.0);
    red_distr_.resize(grid.dimension() * module.numDirection(), 0.0);
    blue_distr_.resize(grid.dimension() * module.numDirection(), 0.0);
//    distr_.resize(grid.dimension(),std::vector<double>(module.numDirection(), 0.0));
}
