#ifndef INITSTATE_HEADER_INCLUDED
#define INITSTATE_HEADER_INCLUDED
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>
#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/utility/SimulatorState.hpp>

void
initState(SimulatorState& state, 
          const GridManager& grid, 
          const FluidProperties& red, 
          const FluidProperties& blue, 
          const LatticeBoltzmannModule& module);
#endif //INITSTATE_HEADER_INCLUDED
