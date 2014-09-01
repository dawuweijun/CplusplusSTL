#include <opm/lattice/utility/SimulatorState.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>

void SimulatorState::init(const GridManager& grid, const LatticeBoltzmannModule& module)
{
    num_phases_ = 2;
    press_.resize(grid.dimension(), 0.0);
    red_density_.resize(grid.dimension(), 0.0);
    blue_density_.resize(grid.dimension(), 0.0);
    red_dist_.resize(grid.dimension()*module.numDirection(), 0.0);
    blue_dist_.resize(grid.dimension()*module.numDirection(), 0.0);
}
