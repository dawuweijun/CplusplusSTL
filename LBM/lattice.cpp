#include "lattice.hpp"
#include <vector>
Lattice::Lattice(const Grid& grid, const Module& module, const std::vector<Fluid>& fluid)
        :externalForce_(fluid[0].tau()*(grid.spaceDim()+1)*4/(module.numDirection()-1))
        ,fluxForce_(0.0)
        ,flux_(0.0)
        ,gff_(std::vector<double>(module.numDirection()))
        ,gfs_(std::vector<double>(module.numDirection()))
        ,grid_(grid)
        ,module_(module)
        ,fluid_(fluid)
{

}


Lattice::~Lattice()
{

}

void
Lattice::potential()
{
    /* Shan-Chen "potential" strengths */
    /* Care must be taken when these parameters are set! */
    const double G0 = 0.0;
    const double G1 = -0.2;
    const double G2 = -0.1;
    const double GS = -0.01;
    double GFF[19] = {G1, G1, G1, G1, G1, G1, 
            		    G2, G2, G2, G2, G2, G2,
		                G2, G2, G2, G2, G2, G2,G0};
  
    double GFS[19] = {GS, GS, GS, GS, GS, GS, 
            		    GS, GS, GS, GS, GS, GS,
		                GS, GS, GS, GS, GS, GS, 
        	    	    G0};
    for (int i = 0; i < gff_.size(); ++i) {
        gff_[i] = GFF[i];
        gfs_[i] = GFS[i];     
    }
}

void
Lattice::propagationBySwap(const Grid& grid, Fluid& fluid)
{
    //swap locally
    for (int i = 0; i <: )
}
void
Lattice::collision(const Grid& grid, const std::vector<Fluid>& fluid)
{

}

void
Lattice::streamingSwap(const Grid& grid, const std::vector<Fluid>& fluid)
{

}
