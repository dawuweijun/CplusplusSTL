///Program for initilization with capillary pressure for SPE9 contest.
///
///$ID: miliu, Thu Aug 15 16:53:36 CST 2013 exp$
///$Email:  miliu@statoil.com$
///

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H
#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorCompressibleTwophase.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <boost/scoped_ptr.hpp>
#include <opm/sample/utils/CapillaryBalance.hpp> 
namespace
{
    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }
} // anon namespace

//--------------------------Main program------------------------
/// Main function initializa the spe9 test
int 
main (int argc, char* argv[])
{
	using namespace Opm;

	std::cout << "\n============== Test program for initilizing spe9 test ===========\n\n";
	std::cout << "Usage:\n";
	std::cout << "  deck_filename=yourfile  // this is your input file\n";
	parameter::ParameterGroup param(argc, argv, false);
	// create the deck_filename to read
	bool use_deck = param.has("deck_filename");
	boost::scoped_ptr<EclipseGridParser> deck;
	boost::scoped_ptr<GridManager> grid;	
	if (use_deck) {
		std::string deck_filename = param.get<std::string>("deck_filename");
		std::string output = "output" + deck_filename;
		deck.reset(new EclipseGridParser(deck_filename));
		grid.reset(new GridManager(*deck));
		//
		Opm::CapillaryBalance cap;
		
		std::ifstream is;
		is.close();
		is.clear();
		is.open(deck_filename.c_str());
		std::ofstream os(output.c_str());
		
		std::cout << "Start initlize the state for SPE9" << std::endl;
		std::cout << "Number of cells: " << grid->c_grid()->number_of_cells << std::endl;
		cap.initState(*grid->c_grid(), is);
		std::cout << "End of initState" << std::endl;
		 
		cap.output(os);
		//output vtkfile
		std::ofstream vtkfile("initState.vtu");
		DataMap dm;
		dm["pressure"] = &cap.pressure();
		dm["swat"] = &cap.swat();
		dm["soil"] = &cap.soil();
		dm["sgas"] = &cap.sgas();
		Opm::writeVtkData(*grid->c_grid(), dm, vtkfile);
	} else {
        warnIfUnusedParams(param);
	}

	return 0;
}
