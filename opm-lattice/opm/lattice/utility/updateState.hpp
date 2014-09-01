#ifndef UPDATESTATE_HEADER_INCLUDED
#define UPDATESTATE_HEADER_INCLUDED
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/utility/SimulatorState.hpp>
#include <opm/lattice/FluidProperties.hpp>
void
updateState(SimulatorState& state, 
            const GridManager& grid, 
            const FluidProperties& red, 
            const FluidProperties& blue, 
            const LatticeBoltzmannModule& module);
#endif //UPDATESTATE_HEADER_INCLUDED
