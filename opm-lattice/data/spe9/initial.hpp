//This program for initializing the pressure
//and saturation distribution for SPE 9 contest
//
//
//$ID: miliu Fri Aug  9 09:47:41 CST 2013 exp$
//$Email: miliu@statoil.com$
//
//
#ifndef OPM_INITIAL_HEADER_INCLUDE
#define OPM_INITIAL_HEADER_INCLUDE

struct UnstructuredGrid;

namespace Opm
{
	namespace parameter { class  ParameterGroup; }
	class EclipseGridParser;
	class IncomPropertiesInterface;
	class BlackoilPropertiesInterface;
	class TwoPhaseState;

	class Initial {
		public:
			Initial(double woc_depth){};
			void read(std::ifstream is,
					  std::vector <double> swat,
					  std::vector <double> pcwo);
			// Find the internal number which the value
			// just placed in.
			//
			int findNu(const double in
					   const std::vector<double> swat,
					   int num);

			double interpolation(const int num,
								 const std::vector<double> pcwo,
								 double result);	
			//
			//
			// initialize oil phase pressure by gravity
			std::vector <double>
			initialPo(const UnstruturedGrid& grid,
					  const BlackoilPropertiesInterface& props,
				      const parameter::ParameterGroup& param,
					  const double gravity);


			// initialize water phase pressure by caillary and gravity
			std::vector <double>
			initialPw(const UnstructuredGrid& grid,
					  const BlackoilPropertiesInterface& props,
					  const parameter::ParameterGroup& param,
					  const double gravity);


			// initialize water saturation by capillary and reserve oil
			// distribution
			std::vecrot <double>
			initialSw(const UnstructuredGrid& grid,
					  const BlackoilPropertiesInterface& props,
					  const parameter::ParameterGroup& param);
			// output the initial distribution to the file 
			// by keyword "PRESSURE SWAT SOIL SGAS "
			void printfvalues(const UnstructuredGrid& grid,
							 const BlackoilPropertiesInterface& props
							 std::ofstream os);

			~Initial();
		private:

			double woc_depth;
			double oil_res;
			
			std::vecror <double> soil;
			std::vector <double> sgas;
			std::vector <double> swat;
			std::vector <double> pcwo;
						
	}
}
