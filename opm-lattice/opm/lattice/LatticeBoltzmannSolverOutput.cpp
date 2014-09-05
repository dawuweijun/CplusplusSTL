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
    std::ostringstream vtkbluename, vtkredname, vtkpname;
    vtkbluename << output_dir << "/vtk_files";
    vtkredname << output_dir << "/vtk_files";
    vtkpname << output_dir << "/vtk_files";
    boost::filesystem::path fpath(vtkbluename.str());
    try {
        create_directories(fpath);
    }
    catch (...) {
        std::cout << "Creating directories failed: " << fpath << std::endl;
        exit(1);
    }
    vtkbluename << "/blue-" << std::setw(6) << std::setfill('0') << step << ".vti";
    std::ofstream vtkblue(vtkbluename.str().c_str());
    vtkredname << "/red-" << std::setw(6) << std::setfill('0') << step << ".vti";
    std::ofstream vtkred(vtkredname.str().c_str());
    vtkpname << "/press-" << std::setw(6) << std::setfill('0') << step << ".vti";
    std::ofstream vtkpress(vtkpname.str().c_str());
    if (!vtkblue && !vtkred) {
        std::cout << "Failed to open " << vtkbluename.str() << vtkredname.str() << vtkpname.str()<<std::endl;
        exit(1);
    }
    DataMap dm, ddm, dmp;
    dmp["pressure"] = &state.pressure();
    dm["reddensity"] = &state.redDensity();
    ddm["bluedensity"] = &state.blueDensity();
    writeVtkData(grid, dm, vtkred);
    writeVtkData(grid, ddm, vtkblue);
    writeVtkData(grid, dmp, vtkpress);
}


void 
outputStateMatlab(const GridManager& grid,
                  const SimulatorState& state,
                  const int step,
                  const std::string& output_dir)
{
    DataMap dm;
    dm["reddensity"] = &state.redDensity();
    dm["bluedensity"] = &state.blueDensity();
    dm["pressure"] = &state.pressure();
    dm["planepressure"] = &state.planePressure();
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
//        fname << "/" << std::setw(6) << std::setfill('0') << step << ".txt";
        fname << "/" << step << ".txt";
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
