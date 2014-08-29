#include <opm/lattice/utility/SimulatorState.hpp>
#include <opm/lattice/GridManager.hpp>

void SimulatorState::init(const GridManager& grid)
{
    num_phases_ = 2;
    press_.resize(grid.dimension(), 0.0);
    vel_.resize(num_phases_ * grid.dimension(), 0.0);
    red_density_.resize(grid.dimension(), 0.0);
    blue_density_.resize(grid.dimension(), 0.0);
}
