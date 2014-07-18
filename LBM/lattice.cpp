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
    for (int i = 0; i < fluid.distribution().size(); ++i) {
        for (int k = 0; k < fluid.distribution()[i].size()-1;) {
            std::swap(fluid.distribution()[i][k], fluid.distribution()[i][k+1]);
            k = k + 2;
        }
    }
    const int NX = grid.NX();
    const int NY = grid.NY();
    const int NZ = grid.NZ();
    //swap between neigbours
    for(int z = 0; z < NZ; ++z) {
		double ztop = (z >= NZ-1 ? z+1-NZ : z+1);
	    double zbot = (z <= 0 ? z-1+NZ : z-1);
        for(int y = 0; y < NY; ++y) {
	        double ytop = (y >= NY-1 ? y+1-NY : y+1);
        	double ybot = (y <= 0 ? y-1+NY : y-1);
            for(x = 0; x < NX; ++x) {
	            xtop = (x >= NX-1 ? x+1-NX : x+1);
            	xbot = (x <= 0 ? x-1+NX : x-1);
                int idx1 = grid.index(x, y, z);
                int idx2 = grid.index(xbot, y, z);
                std::swap(fluid.distribution[idx1][0], fluid.distribution[idx2][1]);
                idx2 = grid.index(x, ybot, z);
                std::swap(fluid.distribution[idx1][2], fluid.distribution[idx2][3]);
                idx2 = grid.index(x, y, zbot);
                std::swap(fluid.distribution[idx1][4], fluid.distribution[idx2][5]);
                idx2 = grid.index(xbot, ybot, z);
                std::swap(fluid.distribution[idx1][6], fluid.distribution[idx2][7]);
                std::swap(fluid.distribution[idx1][8], fluid.distribution[idx2][9]);
                idx2 = grid.index(xbot, y, zbot);
                std::swap(fluid.distribution[idx1][10], fluid.distribution[idx2][11]);
                std::swap(fluid.distribution[idx1][12], fluid.distribution[idx2][13]);
                idx2 = grid.index(x, ybot, zbot);
                std::swap(fluid.distribution[idx1][14], fluid.distribution[idx2][15]);
                std::swap(fluid.distribution[idx1][16], fluid.distribution[idx2][17]);
	        }
        }
    }
}
void
Lattice::collision(const Grid& grid, const std::vector<Fluid>& fluid)
{

}

void
Lattice::streamingSwap(const Grid& grid, const std::vector<Fluid>& fluid)
{

}
