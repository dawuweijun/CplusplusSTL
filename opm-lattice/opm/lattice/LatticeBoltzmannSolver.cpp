#include <opm/lattice/LatticeBoltzmannSolver.hpp>
#include <opm/lattice/LatticeBoltzmannSolverOutput.hpp>
#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/GridManager.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>

#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <numeric>
#include <iomanip>
#include <algorithm>

LatticeBoltzmannSolver::LatticeBoltzmannSolver(const GridManager& grid, const LatticeBoltzmannModule& module, const FluidProperties& red, const FluidProperties& blue)
        :externalForce_(red.tau()*(grid.spaceDim()+1)/((module.numDirection()-1)*2.0))
        ,fluxForce_(0.0)
        ,grid_(grid)
        ,module_(module)
        ,red_(red)
        ,blue_(blue)
{
}

LatticeBoltzmannSolver::SolutionState::SolutionState(const int size)
    : red(std::vector<double>(size))
    , blue(std::vector<double>(size))
{
}

void
LatticeBoltzmannSolver::step(const double dt, SimulatorState& x)
{
    const int size = grid_.dimension()*module_.numDirection();
    SolutionState state(size); 
    state = initVariables(x);
    collisionStepScRed(state);
    collisionStepScBlue(state);
    streamingSwap(state);
    massMomentumCalc(state, x);
    pressureCalc(state, x);
    updateState(state, x);
}


LatticeBoltzmannSolver::SolutionState
LatticeBoltzmannSolver::initVariables(SimulatorState& x)
{
    const int N = grid_.dimension();
    const int ND = module_.numDirection();
    SolutionState state(N*ND);
    
    assert(x.redDist().size() == state.red.size());
    assert(x.blueDist().size() == state.blue.size());

    std::copy(x.redDist().begin(), x.redDist().end(), state.red.begin()); 
    std::copy(x.blueDist().begin(), x.blueDist().end(), state.blue.begin()); 

    return state;
}

void
LatticeBoltzmannSolver::updateState(const SolutionState& state, SimulatorState& x)
{
    // update density.
    assert(x.redDist().size() == state.red.size());
    assert(x.blueDist().size() == state.blue.size());

    std::copy(state.red.begin(), state.red.end(), x.redDist().begin()); 
    std::copy(state.blue.begin(), state.blue.end(), x.blueDist().begin()); 
}

void
LatticeBoltzmannSolver::propagationBySwap(std::vector<double>& dis)
{
    //swap locally
    const int NX = grid_.NX();
    const int NY = grid_.NY();
    const int NZ = grid_.NZ();
    const int N = grid_.dimension();
    const int ND = module_.numDirection();
    assert(dis.size() == N*ND);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < ND-1; ) {
            std::swap(dis[i*ND + k], dis[i*ND + (k+1)]);
            k = k + 2;
        }
    }
    //swap between neigbours
    for(int z = 0; z < NZ; ++z) {
		const int ztop = (z >= NZ-1 ? z+1-NZ : z+1);
	    const int zbot = (z <= 0 ? z-1+NZ : z-1);
        for(int y = 0; y < NY; ++y) {
	        const int ytop = (y >= NY-1 ? y+1-NY : y+1);
        	const int ybot = (y <= 0 ? y-1+NY : y-1);
            for(int x = 0; x < NX; ++x) {
//	            int xtop = (x >= NX-1 ? x+1-NX : x+1);
            	const int xbot = (x <= 0 ? x-1+NX : x-1);
                const int idx1 = grid_.index(x, y, z);
                int idx2 = grid_.index(xbot, y, z);
                std::swap(dis[idx1*ND + 0], dis[idx2*ND + 1]);
                idx2 = grid_.index(x, ybot, z);
                std::swap(dis[idx1*ND + 2], dis[idx2*ND + 3]);
                idx2 = grid_.index(x, y, zbot);
                std::swap(dis[idx1*ND + 4], dis[idx2*ND + 5]);
                idx2 = grid_.index(xbot, ybot, z);
                std::swap(dis[idx1*ND + 6], dis[idx2*ND + 7]);
                idx2 = grid_.index(xbot, ytop, z);
                std::swap(dis[idx1*ND + 8], dis[idx2*ND + 9]);
                idx2 = grid_.index(xbot, y, zbot);
                std::swap(dis[idx1*ND + 10], dis[idx2*ND + 11]);
                idx2 = grid_.index(xbot, y, ztop);
                std::swap(dis[idx1*ND + 12], dis[idx2*ND + 13]);
                idx2 = grid_.index(x, ybot, zbot);
                std::swap(dis[idx1*ND + 14], dis[idx2*ND + 15]);
                idx2 = grid_.index(x, ybot, ztop);
                std::swap(dis[idx1*ND + 16], dis[idx2*ND + 17]);
	        }
        }
    }
}

void
LatticeBoltzmannSolver::collisionStepScBlue(SolutionState& state)
{
    const int N = grid_.dimension();
    const int ND = module_.numDirection();
    double RedDensity, BlueDensity;
    std::vector<double> f(N*ND, 0.0);
    const std::vector<int>& solid = grid_.boundary();
    const std::vector<double>& cx = module_.xVelocity();
    const std::vector<double>& cy = module_.yVelocity();
    const std::vector<double>& cz = module_.zVelocity();
    const std::vector<double>& w = module_.weight();
    const std::vector<double>& gff = blue_.gff();
    const std::vector<double>& gfs = blue_.gfs();
    for (int i = 0; i < N; ++i) {
        if (solid[i] == 0) {
            RedDensity = 0.0;
            for (int k = 0; k < ND; ++k) {
                RedDensity += state.red[i*ND + k];
            }
            for (int k = 0; k < ND; ++k) {
                f[i*ND + k] = RedDensity;
            }
        } else {
            for (int k = 0; k < ND; ++k) {
                f[i*ND + k] = -10.0;
            }
        }
    }
    propagationBySwap(f);
    fcalcSc(f, gff, gfs);
    //
    //Bounce-back
    //
    std::vector<double> RedAfter(ND, 0.0);
    std::vector<double> BlueAfter(ND, 0.0);
    for (int i = 0; i < N; ++i) {
        if (solid[i] == 1) {
            for (int k = 0; k < ND-1; ++k) {
                if (k%2 == 0) {
                    BlueAfter[k] = state.blue[i*ND + (k+1)];
                } else {
                    BlueAfter[k] = state.blue[i*ND + (k-1)];
                }
            }
            for (int k = 0; k < ND-1; ++k) {
                state.blue[i*ND + k] = BlueAfter[k];
            }
        } else {
            std::vector<double> velocity(grid_.spaceDim(), 0.0);
            RedDensity = BlueDensity=0;
            for (int k = 0; k < ND; ++k) {
                RedDensity += state.red[i*ND + k];
                BlueDensity += state.blue[i*ND + k];
                velocity[0] += (state.red[i*ND + k] + state.blue[i*ND + k])*cx[k] 
                               - (1.0 / blue_.lambda())*(std::fabs(state.blue[i*ND + k]*fluxForce_*cx[k]));
                velocity[1] += (state.red[i*ND + k] + state.blue[i*ND + k])*cy[k];
                velocity[2] += (state.red[i*ND + k] + state.blue[i*ND + k])*cz[k];
            }
            //add Sc potential to momentum.
            for (int k = 0; k < 3; ++k) {
                velocity[k] += -(1.0 / blue_.lambda()) * BlueDensity * f[i*ND + k];
            }
            for (int k = 0; k < ND; ++k) {
                state.blue[i*ND + k] += blue_.lambda()*(state.blue[i*ND + k] - NipSc(0, RedDensity, BlueDensity, cx[k], cy[k], cz[k], w[k], velocity));
            }
        }
    }
}

void
LatticeBoltzmannSolver::collisionStepScRed(SolutionState& state)
{
    const int N = grid_.dimension();
    const int ND = module_.numDirection();
    double RedDensity, BlueDensity;
    std::vector<double> f(N*ND, 0.0);
    const std::vector<int>& solid = grid_.boundary();
    const std::vector<double>& cx = module_.xVelocity();
    const std::vector<double>& cy = module_.yVelocity();
    const std::vector<double>& cz = module_.zVelocity();
    const std::vector<double>& w = module_.weight();
    const std::vector<double>& gff = red_.gff();
    const std::vector<double>& gfs = red_.gfs();
    for (int i = 0; i < N; ++i) {
        if (solid[i] == 0) {
            BlueDensity = 0.0;
            for (int k = 0; k < ND; ++k) {
                BlueDensity += state.blue[i*ND + k];
            }
            for (int k = 0; k < ND; ++k) {
                f[i*ND + k] = BlueDensity;
            }
        } else {
            for (int k = 0; k < ND; ++k) {
                f[i*ND + k] = -10.0;
            }
        }
    }
    propagationBySwap(f);
    fcalcSc(f, gff, gfs);
    //
    //Bounce-back
    //
    std::vector<double> RedAfter(ND, 0.0);
    std::vector<double> BlueAfter(ND, 0.0);
    for (int i = 0; i < N; ++i) {
        if (solid[i] == 1) {
            for (int k = 0; k < ND-1; ++k) {
                if (k%2 == 0) {
                    RedAfter[k] = state.red[i*ND + (k+1)];
                } else {
                    RedAfter[k] = state.red[i*ND + (k-1)];
                }
            }
            for (int k = 0; k < ND-1; ++k) {
                state.red[i*ND + k] = RedAfter[k];
            }
        } else {
            std::vector<double> velocity(grid_.spaceDim());
            RedDensity = BlueDensity=0;
            for (int k = 0; k < ND; ++k) {
                RedDensity += state.red[i*ND + k];
                BlueDensity += state.blue[i*ND + k];
                velocity[0] += (state.red[i*ND + k] + state.blue[i*ND + k])*cx[k] 
                               - (1.0 / red_.lambda())*(std::fabs(state.red[i*ND + k]*fluxForce_*cx[k]));
                velocity[1] += (state.red[i*ND + k] + state.blue[i*ND + k])*cy[k];
                velocity[2] += (state.red[i*ND + k] + state.blue[i*ND + k])*cz[k];
            }
            //add Sc potential to momentum.
            for (int k = 0; k < 3; ++k) {
                velocity[k] += -(1.0 / red_.lambda())*RedDensity*f[i*ND + k];
            }

            for (int k = 0; k < ND; ++k) {
                state.red[i*ND + k] += red_.lambda()*(state.red[i*ND + k] - NipSc(1, RedDensity, BlueDensity, cx[k], cy[k], cz[k], w[k], velocity));
            }
        }
    }
}

double
LatticeBoltzmannSolver::NipSc(const int flag, const double Rden, const double Bden, const double cxk, const double cyk, const double czk, const double wk, std::vector<double> velocity)
{
    double Nip = 0, uc = 0, u2 = 0;
    if ((Rden+Bden) > 0.0) {
        for (int i = 0; i < 3; ++i) {
            velocity[i] /= Rden+Bden;
        }
    }
    uc = velocity[0]*cxk + velocity[1]*cyk + velocity[2]*czk;
    for (int i = 0; i < 3; ++i) {
        u2 += std::pow(velocity[i], 2);
    }
    Nip = wk*(1 + 3*uc + 4.5*uc*uc - 1.5*u2);
    return flag==0?Bden*Nip:Rden*Nip;
}

void
LatticeBoltzmannSolver::fcalcSc(std::vector<double>& f, 
                                const std::vector<double>& gff, 
                                const std::vector<double>& gfs)
{
    const int N = grid_.dimension();
    const int ND = module_.numDirection();
    assert(f.size() == N*ND);
    const std::vector<double>& cx = module_.xVelocity();
    const std::vector<double>& cy = module_.yVelocity();
    const std::vector<double>& cz = module_.zVelocity();
    for (int i = 0; i < N; ++i) {
        double fx = 0, fy = 0, fz = 0;
        for (int k = 0; k < ND; ++k) {
            if (f[i*ND + k] > 0.0) {
                fx += gff[k]*cx[k]*f[i*ND + k];
                fy += gff[k]*cy[k]*f[i*ND + k];
                fz += gff[k]*cz[k]*f[i*ND + k];
            } else if (f[i*ND + k] < -9.9) {
                fx += gfs[k]*cx[k];
                fy += gfs[k]*cy[k];
                fz += gfs[k]*cz[k];
            }
        }
        f[i*ND + 0] = -fx;
        f[i*ND + 1] = -fy;
        f[i*ND + 2] = -fz;
    }
}
void
LatticeBoltzmannSolver::streamingSwap(SolutionState& state)
{
    propagationBySwap(state.red);
    propagationBySwap(state.blue);
    
}

void
LatticeBoltzmannSolver::massMomentumCalc(const SolutionState& state, SimulatorState& x)
{
    const int N = grid_.dimension();
    const int ND = module_.numDirection();
    std::vector<double> mass(2, 0.0);
    std::vector<double> press(N, 0.0);
    double maxabsvelocity = -10.0;
    mass[0] = std::accumulate(state.red.begin(), state.red.end(), 0.0);
    mass[1] = std::accumulate(state.blue.begin(), state.blue.end(), 0.0);
    const std::vector<int>& boundary = grid_.boundary();
    const std::vector<double>& cx = module_.xVelocity();
    const std::vector<double>& cy = module_.yVelocity();
    const std::vector<double>& cz = module_.zVelocity();
    for (int i = 0; i < N; ++i) {
        std::vector<double> velocity(3, 0.0);
        double tmp = 0;
        if (boundary[i] == 0) {
            for (int k = 0; k < ND; ++k) {
                tmp += state.red[i*ND + k] + state.blue[i*ND + k];
                velocity[0] += (state.red[i*ND + k] + state.blue[i*ND + k])*cx[k];    
                velocity[1] += (state.red[i*ND + k] + state.blue[i*ND + k])*cy[k];    
                velocity[2] += (state.red[i*ND + k] + state.blue[i*ND + k])*cz[k];    
            }
            for (int k = 0; k < 3; ++k) {
                velocity[k] /= tmp; 
            }
            double velsquare = velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2];
            press[i] = tmp / 3.0*(1.0 - velsquare);
            if (std::sqrt(velsquare) > maxabsvelocity) {
                maxabsvelocity = std::sqrt(velsquare);
            }
        }
    }
    assert(press.size() == x.pressure().size());
    std::copy(press.begin(), press.end(), x.pressure().begin());
    std::cout << "\n   Mass Momentum Summary:\n";
    std::cout << std::setw(18) << std::setprecision(9) << "     Total Mass   : " << mass[0] + mass[1] << std::endl;
    std::cout << std::setw(18) << std::setprecision(9) << "     Red   Mass   : " << mass[0] << std::endl;
    std::cout << std::setw(18) << std::setprecision(9) << "     Blue  Mass   : " << mass[1] << std::endl;
    std::cout << std::setw(18) << std::setprecision(9) << "     FluxForce    : " << fluxForce_ << std::endl;
    std::cout << std::setw(18) << std::setprecision(9) << "     ExternalForce: " << externalForce_ << std::endl;
    std::cout << std::setw(18) << std::setprecision(9) << "     Max abs vel  : " << maxabsvelocity << std::endl;
}

void
LatticeBoltzmannSolver::pressureCalc(const SolutionState& state, SimulatorState& x)
{
    //get plane pressure
    std::vector<double> planepress;
    typedef std::vector<double>::size_type vec_sz;
    for (int xx = 0; xx < grid_.NX(); ++xx) {
        std::vector<int> idx;
        planepress.push_back(xx);
        double p = 0, c = 0;
        for (int y = 0; y < grid_.NY(); ++y) {
            for (int z = 0; z < grid_.NZ(); ++z) {
                idx.push_back(grid_.index(xx, y, z));
            }
        }
        for (vec_sz i = 0; i < idx.size(); ++i) {
            if (state.red[idx[i]] != 0 && grid_.boundary()[idx[i]] != 1) {
                p += x.pressure()[idx[i]];
                c++;
            }
        }
        planepress.push_back(c>0?p/c:p);
    }
    x.planePressure().resize(planepress.size());
    std::copy(planepress.begin(), planepress.end(), x.planePressure().begin());
}
