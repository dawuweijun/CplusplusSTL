#ifndef LATTICEBOLTZMANNSOLVEROUTPUT_INCLUDED_HEADER
#define LATTICEBOLTZMANNSOLVEROUTPUT_INCLUDED_HEADER
#include <opm/lattice/utility/DataMap.hpp>
#include <opm/lattice/utility/writeVtkData.hpp>

void outputStateVtk(const GridManager& grid,
                    const SimulatorState& state,
                    const int step,
                    const std::string& output_dir);
void outputStateMatlab(const GridManager& grid,
                       const SimulatorState& state,
                       const int step,
                       const std::string& output_dir);

#endif // LATTICEBOLTZMANNSOLVEROUTPUT_INCLUDED_HEADER
