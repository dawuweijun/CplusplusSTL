#ifndef WRITEVTKDATA_HEADER_INCLUDED
#define WRITEVTKDATA_HEADER_INCLUDED

#include <opm/lattice/utility/DataMap.hpp>
#include <opm/lattice/GridManager.hpp>
#include <iosfwd>
/// Vtk output for cartesian grids.
void writeVtkData(const GridManager& grid,
                  const DataMap& data,
                  std::ostream& os);

#endif // WRITEVTKDATA_HEADER_INCLUDED
