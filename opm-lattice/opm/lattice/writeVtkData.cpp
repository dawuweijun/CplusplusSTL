#include <opm/lattice/DataMap.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/writeVtkData.hpp>
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
    os.precision(12);
    os << "<?xml version=\"1.0\"?>\n";
    PMap pm;
    pm["type"] = "ImageData";
    Tag vtkfiletag("VTKFile", pm, os);
    Tag ugtag("ImageData", os);
    int num_pts = grid.dimension();
    pm.clear();
    pm["NumberOfPoints"] = boost::lexical_cast<std::string>(num_pts);
    Tag piecetag("Piece", pm, os);
    {
        Tag pointstag("Points", os);
        pm.clear();
        pm["type"] = "Float64";
        pm["Name"] = "Coordinates";
        pm["NumberOfComponents"] = "3";
        pm["format"] = "ascii";
        Tag datag("DataArray", pm, os);
        for (int z  = 0; z < grid.NZ(); ++z) {
            for (int y = 0; y < grid.NY(); ++y) {
                for (int x = 0; x < grid.NX(); ++x) {
                    Tag::indent(os);
                    os << x << ' ' << y << ' ' << z << '\n';
                }
            }
        }
    }
    {
        pm.clear();
//        if (data.find("density") != data.end()) {
//            pm["Scalars"] = "density";
//        } else if (data.find("pressure") != data.end()) {
//            pm["Scalars"] = "pressure";
//        }
        pm["Scalars"] = "density";
        Tag celldatatag("CellData", pm, os);
        pm.clear();
        pm["NumberOfComponents"] = "1";
        pm["format"] = "ascii";
        pm["type"] = "Float64";
        for (DataMap::const_iterator dit = data.begin(); dit != data.end(); ++dit) {
            pm["Name"] = dit->first;
            const std::vector<double>& field = *(dit->second);
            const int num_comps = field.size() / num_pts;
            pm["NumberOfComponents"] = boost::lexical_cast<std::string>(num_comps);
            Tag ptag("DataArray", pm, os);
            const int num_per_line =  num_comps == 1 ? 5 : num_comps;
            for (int item = 0; item < num_pts*num_comps; ++item) {
                if (item % num_per_line == 0) {
                    Tag::indent(os);
                }
                double value = field[item];
                if (std::fabs(value) < std::numeric_limits<double>::min()) {
                    // Avoiding denormal numbers to work around
                    // bug in Paraview.
                    value = 0.0;
                }
                os << value << ' ';
                if (item % num_per_line == num_per_line - 1
                    || item == num_pts - 1) {
                    os << '\n';
                }
            }
        }
    }
}


