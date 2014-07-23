#ifndef WRITEVTKDATA_HEADER_INCLUDED
#define WRITEVTKDATA_HEADER_INCLUDED

#include <string>
#include <map>
#include <vector>
#include <array>
#include <iosfwd>
#include <opm/lattice/DataMap.hpp>
/// Vtk output for cartesian grids.
void writeVtkData(const std::array<int, 3>& dims,
                  const std::array<double, 3>& cell_size,
                  const DataMap& data,
                  std::ostream& os);

#endif // WRITEVTKDATA_HEADER_INCLUDED
