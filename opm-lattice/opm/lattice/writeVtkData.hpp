#ifndef WRITEVTKDATA_HEADER_INCLUDED
#define WRITEVTKDATA_HEADER_INCLUDED

#include <opm/lattice/DataMap.hpp>
#include <opm/lattice/GridManager.hpp>
#include <string>
#include <map>
#include <vector>
#include <array>
#include <iosfwd>
/// Vtk output for cartesian grids.
void writeVtkData(const GridManager& grid,
                  const DataMap& data,
                  std::ostream& os);

#endif // WRITEVTKDATA_HEADER_INCLUDED
