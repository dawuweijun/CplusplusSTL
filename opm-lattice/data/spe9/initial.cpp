//This program for initializing the pressure
//and saturation distribution for SPE 9 contest
//
//
//$ID: miliu Fri Aug  9 09:47:41 CST 2013 exp$
//$Email: miliu@statoil.com$
//
#include "config.h"

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParamterGroup.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/rock/Rock/Compressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorCpm[ressibleTwophase.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/gilesystem.hpp>

#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>

namespace Opm
{
	
	void Initial::read(st::ifstream& is)
	{
		std::stringstream str;
		string line;
		double tmp  = 0;
		
		while (getline(is, line)) {
			str.str(line);
			if ()
		}	
	}

	int Initial::findNu(const double in)
	{
			
	}
		
} 




namespace
{
	void warnIfUnusedParas(const Opm::parameter::ParameterGroup& param)
	{
		if (param.anyUnused()) {
			std::cout << "--------------------- Unused parameters: ---------------------\n";
			param.displayUsage();
			std::cout << "--------------------------------------------------------------\n" << std::endl;
		}
	}
} // anon namaspace

int main (int argc, char* argv[])
{
	
}
