#include <opm/lattice/utility/DataMap.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/utility/writeVtkData.hpp>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <array>
#include <iosfwd>
#include <cmath>
#include <boost/lexical_cast.hpp>

void writeVtkData(const GridManager& grid,
                  const DataMap& data,
                  std::ostream& os)
{
    if (grid.spaceDim() != 3) {
        std::cout << "Vtk output for 3d grids only" << std::endl;
        exit(1);
    }
    const int num_pts = grid.dimension();
    os.precision(12);
    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "  <ImageData WholeExtent=\"0 "<<grid.NX()-1<<" 0 "<< grid.NY()-1<<" 0 "<<grid.NZ()-1<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    os << "  <Piece Extent=\"0 "<<grid.NX()-1<<" 0 "<<grid.NY()-1<<" 0 "<<grid.NZ()-1<<"\">\n";
    os << "    <PointData Scalars=\"density\">\n";
    for (DataMap::const_iterator dit = data.begin(); dit != data.end(); ++dit) {
        const std::vector<double>& field = *(dit->second);
        const int num_comps = field.size() / num_pts;
        os << "       <DataArray type=\"Float32\" Name=\""<<dit->first<<"\" NumberOfComponents=\""<<boost::lexical_cast<std::string>(num_comps)<<"\" format=\"ascii\">\n"; 
        const int num_per_line =  num_comps == 1 ? 5 : num_comps;
        for (int i= 0; i < num_pts; ++i) {
            if (i % num_per_line == 0) {
                os <<"      ";
            }
            double value = field[i];
            if (std::fabs(value) < std::numeric_limits<double>::min()) {
                // Avoiding denormal numbers to work around
                // bug in Paraview.
                value = 0.0;
            }
            os << value << "   ";
            if (i % num_per_line == num_per_line - 1
                || i == num_pts - 1) {
                os << '\n';
            }
        }
        os <<"       </DataArray>\n";
    }
    os << "    </PointData>\n";
    os << "    <CellData>\n";
    os << "    </CellData>\n";
    os << "  </Piece>\n";
    os << " </ImageData>\n";
    os << "</VTKFile>\n";
}
