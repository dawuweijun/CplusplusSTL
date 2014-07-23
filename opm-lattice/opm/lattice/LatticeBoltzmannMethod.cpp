#include <opm/lattice/lattice.hpp>
#include <vector>
LatticeBoltzmannMethod::LatticeBoltzmannMethod(const GridManager& grid, const LatticeBoltzmannModule& module, const FluidProperties& red, const FluidProperties& blue)
        :externalForce_(red.tau()*(grid.spaceDim()+1)*4/(module.numDirection()-1))
        ,fluxForce_(0.0)
        ,flux_(0.0)
        ,gff_(std::vector<double>(module.numDirection()))
        ,gfs_(std::vector<double>(module.numDirection()))
        ,grid_(grid)
        ,module_(module)
        ,red_(red)
        ,blue_(blue)
{

}


LatticeBoltzmannMethod::~LatticeBoltzmannMethod()
{

}

void
LatticeBoltzmannMethod::potential()
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
    for (int i = 0; i < static_cast<int>(gff_.size()); ++i) {
        gff_[i] = GFF[i];
        gfs_[i] = GFS[i];     
    }
}

void
LatticeBoltzmannMethod::propagationBySwap(std::vector<std::vector<double>>& dis)
{
    //swap locally
    for (int i = 0; i < static_cast<int>(dis.size()); ++i) {
        for (int k = 0; k < static_cast<int>(dis[i].size())-1;) {
            std::swap(dis[i][k], dis[i][k+1]);
            k = k + 2;
        }
    }
    const int NX = grid_.NX();
    const int NY = grid_.NY();
    const int NZ = grid_.NZ();
    //swap between neigbours
    for(int z = 0; z < NZ; ++z) {
		double ztop = (z >= NZ-1 ? z+1-NZ : z+1);
	    double zbot = (z <= 0 ? z-1+NZ : z-1);
        for(int y = 0; y < NY; ++y) {
	        double ytop = (y >= NY-1 ? y+1-NY : y+1);
        	double ybot = (y <= 0 ? y-1+NY : y-1);
            for(int x = 0; x < NX; ++x) {
	            int xtop = (x >= NX-1 ? x+1-NX : x+1);
            	int xbot = (x <= 0 ? x-1+NX : x-1);
                int idx1 = grid_.index(x, y, z);
                int idx2 = grid_.index(xbot, y, z);
                std::swap(dis[idx1][0], dis[idx2][1]);
                idx2 = grid_.index(x, ybot, z);
                std::swap(dis[idx1][2], dis[idx2][3]);
                idx2 = grid_.index(x, y, zbot);
                std::swap(dis[idx1][4], dis[idx2][5]);
                idx2 = grid_.index(xbot, ybot, z);
                std::swap(dis[idx1][6], dis[idx2][7]);
                idx2 = grid_.index(xbot, ytop, z);
                std::swap(dis[idx1][8], dis[idx2][9]);
                idx2 = grid_.index(xbot, y, zbot);
                std::swap(dis[idx1][10], dis[idx2][11]);
                idx2 = grid_.index(xbot, y, ztop);
                std::swap(dis[idx1][12], dis[idx2][13]);
                idx2 = grid_.index(x, ybot, zbot);
                std::swap(dis[idx1][14], dis[idx2][15]);
                idx2 = grid_.index(x, ybot, ztop);
                std::swap(dis[idx1][16], dis[idx2][17]);
	        }
        }
    }
}
void
LatticeBoltzmannMethod::collisionStepScBlue(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist)
{
    const int n = grid_.dimension();
    const int nv = module_.numDirection();
    double RedDensity = 0, BlueDensity = 0;
    std::vector<double> velocity(grid_.spaceDim());
    std::vector<std::vector<double>> f;
    const std::vector<int>& solid = grid_.boundary();
    const std::vector<double>& cx = module_.xVelocity();
    const std::vector<double>& cy = module_.yVelocity();
    const std::vector<double>& cz = module_.zVelocity();
    const std::vector<double>& w = module_.weight();
    for (int i = 0; i < n; ++i) {
        if (solid[i] == 0) {
            for (int k = 0; k < nv; ++k) {
                RedDensity += RedDist[i][k];
                BlueDensity += BlueDist[i][k];
            }
            for (int k = 0; k < nv; ++k) {
                f[i][k] = RedDensity;
            }
        } else if (solid[i] == 1) {
            for (int k = 0; k < nv; ++k) {
                f[i][k] = -10.0;
            }
        }
    }
    propagationBySwap(f);
    fcalcSc(f);
    //
    //Bounce-back
    //
    std::vector<double> RedAfter(nv);
    std::vector<double> BlueAfter(nv);
    for (int i = 0; i < n; ++i) {
        if (solid[i] == 1) {
            for (int k = 0; k < nv-1; ++k) {
                if (k%2 == 0) {
                    BlueAfter[k] = BlueDist[i][k+1];
                } else {
                    BlueAfter[k] = BlueDist[i][k-1];
                }
            }
            for (int k = 0; k < nv-1; ++k) {
                BlueDist[i][k] = BlueAfter[k];
            }
        } else if (solid[i] == 0) {
            RedDensity = BlueDensity=0;
            for (int k = 0; k < nv; ++k) {
                RedDensity += RedDist[i][k];
                BlueDensity += BlueDist[i][k];
                velocity[0] += (RedDist[i][k] + BlueDist[i][k]) * cx[k] 
                               - (1.0 / blue_.lambda()) * (std::fabs(RedDist[i][k] * fluxForce_ * cx[k]));
                velocity[1] += (RedDist[i][k] + BlueDist[i][k]) * cy[k];
                velocity[2] += (RedDist[i][k] + BlueDist[i][k]) * cz[k];
            }
            //add Sc potential to momentum.
            for (int k = 0; k < 3; ++k) {
                velocity[k] += -(1.0 / blue_.lambda()) * BlueDensity * f[i][k];
            }
            for (int k = 0; k < nv; ++k) {
                BlueDist[i][k] += blue_.lambda()*(BlueDist[i][k] - NipSc(0, RedDensity, BlueDensity, cx[k], cy[k], cz[k], w[k], velocity));
            }
        }
    }
}

void
LatticeBoltzmannMethod::collisionStepScRed(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist)
{
    const int n = grid_.dimension();
    const int nv = module_.numDirection();
    double RedDensity = 0, BlueDensity = 0;
    std::vector<double> velocity(grid_.spaceDim());
    std::vector<std::vector<double>> f;
    const std::vector<int>& solid = grid_.boundary();
    const std::vector<double>& cx = module_.xVelocity();
    const std::vector<double>& cy = module_.yVelocity();
    const std::vector<double>& cz = module_.zVelocity();
    const std::vector<double>& w = module_.weight();
    for (int i = 0; i < n; ++i) {
        if (solid[i] == 0) {
            for (int k = 0; k < nv; ++k) {
                RedDensity += RedDist[i][k];
                BlueDensity += BlueDist[i][k];
            }
            for (int k = 0; k < nv; ++k) {
                f[i][k] = RedDensity;
            }
        } else if (solid[i] == 1) {
            for (int k = 0; k < nv; ++k) {
                f[i][k] = -10.0;
            }
        }
    }
    propagationBySwap(f);
    fcalcSc(f);
    //
    //Bounce-back
    //
    std::vector<double> RedAfter(nv);
    std::vector<double> BlueAfter(nv);
    for (int i = 0; i < n; ++i) {
        if (solid[i] == 1) {
            for (int k = 0; k < nv-1; ++k) {
                if (k%2 == 0) {
                    RedAfter[k] = RedDist[i][k+1];
                } else {
                    RedAfter[k] = RedDist[i][k-1];
                }
            }
            for (int k = 0; k < nv-1; ++k) {
                RedDist[i][k] = RedAfter[k];
            }
        } else if (solid[i] == 0) {
            RedDensity = BlueDensity=0;
            for (int k = 0; k < nv; ++k) {
                RedDensity += RedDist[i][k];
                BlueDensity += BlueDist[i][k];
                velocity[0] += (RedDist[i][k] + BlueDist[i][k]) * cx[k] 
                               - (1.0 / red_.lambda()) * (std::fabs(RedDist[i][k] * fluxForce_ * cx[k]));
                velocity[1] += (RedDist[i][k] + BlueDist[i][k]) * cy[k];
                velocity[2] += (RedDist[i][k] + BlueDist[i][k]) * cz[k];
            }
            //add Sc potential to momentum.
            for (int k = 0; k < 3; ++k) {
                velocity[k] += -(1.0 / red_.lambda()) * RedDensity * f[i][k];
            }
            for (int k = 0; k < nv; ++k) {
                BlueDist[i][k] += red_.lambda()*(RedDist[i][k] - NipSc(1, RedDensity, BlueDensity, cx[k], cy[k], cz[k], w[k], velocity));
            }
        }
    }
}

double
LatticeBoltzmannMethod::NipSc(const int flag, const double Rden, const double Bden, const double cxk, const double cyk, const double czk, const double wk, std::vector<double>& velocity)
{
    double Nip = 0, rho = 0, uc = 0, u2 = 0;
    rho = Rden + Bden;
    if (rho <= 0.0) {
        rho = 0;
    } else {
        for (int i = 0; i < static_cast<int>(velocity.size()); ++i) {
            velocity[i] /= rho;
        }
    }
    uc = velocity[0] * cxk + velocity[1] * cyk + velocity[2] * czk;
    for (int i = 0; i < static_cast<int>(velocity.size()); ++i) {
        u2 += std::pow(velocity[i], 2);
    }
    if (flag == 0) {
        Nip = Bden * wk * (1 + 3 * uc + 4.5 * uc * uc - 1.5 * u2); 
    }
    if (flag == 1) {
        Nip = Rden * wk * (1 + 3 * uc + 4.5 * uc * uc - 1.5 * u2); 
    }
    return Nip;
}
void
LatticeBoltzmannMethod::fcalcSc(std::vector<std::vector<double>>& f)
{
    double fx = 0, fy = 0, fz = 0;
    const int n = grid_.dimension();
    const int nv = module_.numDirection();
    const std::vector<double>& cx = module_.xVelocity();
    const std::vector<double>& cy = module_.yVelocity();
    const std::vector<double>& cz = module_.zVelocity();
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < nv; ++k) {
            if (f[i][k] > 0.0) {
                fx += gff_[k] * cx[k] * f[i][k];
                fy += gff_[k] * cy[k] * f[i][k];
                fz += gff_[k] * cz[k] * f[i][k];
            } else if (f[i][k] < -9.9) {
                fx += gfs_[k] * cx[k];
                fy += gfs_[k] * cy[k];
                fz += gfs_[k] * cz[k];
                
            }
        }
        f[i][0] = -fx;
        f[i][1] = -fy;
        f[i][2] = -fz;
    }
}
void
LatticeBoltzmannMethod::streamingSwap(std::vector<std::vector<double>>& RedDist, std::vector<std::vector<double>>& BlueDist)
{
    propagationBySwap(RedDist);
    propagationBySwap(BlueDist);
    
}
