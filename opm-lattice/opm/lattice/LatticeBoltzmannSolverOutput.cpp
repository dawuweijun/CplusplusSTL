#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/utility/DataMap.hpp>
#include <opm/lattice/utility/SimulatorState.hpp>
#include <opm/lattice/utility/writeVtkData.hpp>
#include <boost/filesystem.hpp>
#include <numeric>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
void 
outputStateVtk(const GridManager& grid,
               const SimulatorState& state,
               const int step,
               const std::string& output_dir)
{
    std::ostringstream vtkfilename;
    vtkfilename << output_dir << "/vtk_files";
    boost::filesystem::path fpath(vtkfilename.str());
    try {
        create_directories(fpath);
    }
    catch (...) {
        std::cout << "Creating directories failed: " << fpath << std::endl;
        exit(1);
    }
    vtkfilename << "/output-" << std::setw(6) << std::setfill('0') << step << ".vti";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
        std::cout << "Failed to open " << vtkfilename.str() << std::endl;
        exit(1);
    }
    DataMap dm;
    dm["density"] = &state.velocity();
//    dm["pressure"] = &state.pressure();
    writeVtkData(grid, dm, vtkfile);
}


void 
outputStateMatlab(const GridManager& grid,
                  const SimulatorState& state,
                  const int step,
                  const std::string& output_dir)
{
    DataMap dm;
    dm["density"] = &state.velocity();
    dm["pressure"] = &state.pressure();
    // Write data (not grid) in Matlab format
    for (DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
        std::ostringstream fname;
        fname << output_dir << "/" << it->first;
        boost::filesystem::path fpath = fname.str();
        try {
            create_directories(fpath);
        }
        catch (...) {
            std::cout << "Creating directories failed: " << fpath << std::endl;
            exit(1);
        }
        fname << "/" << std::setw(6) << std::setfill('0') << step << ".txt";
        std::ofstream file(fname.str().c_str());
        if (!file) {
            std::cout << "Failed to open " << fname.str() << std::endl;
            exit(1);
        }
        file.precision(15);
        const std::vector<double>& d = *(it->second);
        std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
    }
}
