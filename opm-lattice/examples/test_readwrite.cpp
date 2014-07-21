//Test Program for reading and writing the data file
//
//
//
//$ID: miliu, Mon Aug 12 15:20:32 CST 2013 exp$
//
//$Email:  miliu@statoil.com$
//


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


#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

#include <opm/sample/utils/initialBasic.hpp> 
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
int 
main (int argc, char* argv[])
{
	using namespace Opm;

	std::cout << "\n============== Test program for reading and writing data file ===========\n\n";
	parameter::ParameterGroup param(argc, argv, false);
	// create the deck_filename to read
	//
	bool use_deck = param.has("deck_filename");
	std::unique_ptr<EclipseGridParser> deck;
		
	if (use_deck) {
		std::string deck_filename = param.get<std::string>("deck_filename");
		std::string output = "output" + deck_filename;
		Opm::InitialBasic init;
		std::ifstream is(deck_filename.c_str());
		std::ofstream os(output.c_str());
		init.classifykeywords(is);
		init.write("SWAT");
		init.write("PCWO");
		init.write("DENSITY");
		init.write("WOC");
		init.write("RESOIL");
		std::cout << "write them to file " << output << std::endl;
		init.write(os);
	} else {
        warnIfUnusedParams(param);
	}
}
