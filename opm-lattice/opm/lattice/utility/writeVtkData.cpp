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


typedef std::map<std::string, std::string> PMap;


struct Tag
{
    Tag(const std::string& tag, const PMap& props, std::ostream& os)
        : name_(tag), os_(os)
    {
        indent(os);
        os << "<" << tag;
        for (PMap::const_iterator it = props.begin(); it != props.end(); ++it) {
            os << " " << it->first << "=\"" << it->second << "\"";
        }
        os << ">\n";
        ++indent_;
    }
    Tag(const std::string& tag, std::ostream& os)
        : name_(tag), os_(os)
    {
        indent(os);
        os << "<" << tag << ">\n";
        ++indent_;
    }
    ~Tag()
    {
        --indent_;
        indent(os_);
        os_ << "</" << name_ << ">\n";
    }
    static void indent(std::ostream& os)
    {
        for (int i = 0; i < indent_; ++i) {
            os << "  ";
        }
    }
private:
    static int indent_;
    std::string name_;
    std::ostream& os_;
};

int Tag::indent_ = 0;


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
    os << "<!-- LBMflow v1.0.1, www.lbmflow.com -->\n";
    os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "  <ImageData WholeExtent=\"0 "<<grid.NX()-1<<" 0 "<< grid.NY()-1<<" 0 "<<grid.NZ()-1<<"\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    os << "  <Piece Extent=\"0 "<<grid.NX()-1<<" 0 "<<grid.NY()-1<<" 0 "<<grid.NZ()-1<<"\">\n";
    os << "    <PointData Scalars=\"pressure\" Vectors=\"density\">\n";
//    os << "       <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"2\" format=\"ascii\">\n"
    for (DataMap::const_iterator dit = data.begin(); dit != data.end(); ++dit) {
        const std::vector<double>& field = *(dit->second);
        const int num_comps = field.size() / num_pts;
        os << "       <DataArray type=\"Float32\" Name=\""<<dit->first<<"\" NumberOfComponents=\""<<boost::lexical_cast<std::string>(num_comps)<<"\" format=\"ascii\">\n"; 
        const int num_per_line =  num_comps == 1 ? 5 : num_comps;
        for (int item = 0; item < num_pts*num_comps; ++item) {
            if (item % num_per_line == 0) {
                os <<"      ";
            }
            double value = field[item];
            if (std::fabs(value) < std::numeric_limits<double>::min()) {
                // Avoiding denormal numbers to work around
                // bug in Paraview.
                value = 0.0;
            }
            os << value << "  ";
            if (item % num_per_line == num_per_line - 1
                || item == num_pts - 1) {
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


