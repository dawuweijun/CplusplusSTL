#include <vector>
#include <iostream>
#include <set>
#include "grid.hpp"

Grid::Grid(double NX, double NY, double NZ)
     :NX_(NX)
     ,NY_(NY)
     ,NZ_(NZ)
     ,boundary_(NX*NY*NZ, 0)
{
    std::cout << "Grid was successfully created!\n";
}

Grid::~Grid()
{}

double
Grid::dimension() const
{
    return NX_ * NY_ * NZ_;
}
double
Grid::spaceDim() const
{
    if (NZ_ == 1)
        return 2;
    else
        return 3;
}
double 
Grid::index(double x, double y, double z) const
{
    return NZ_ * NY_ * z +  NY_ * y + x;
}
std::vector<int> Grid::boundary()
{
    std::vector<int> idx;
    //compute innner index.
    for (int z  = 0; z < NZ_; ++z) {
        for (int y = 0; y < NY_; ++y) {
            idx.push_back(index(0, y, z));
            idx.push_back(index(NX_-1, y, z));
        }
    }
    for (int z = 0; z < NZ_; ++z) {
        for (int x = 0; x < NX_; ++x) {
            idx.push_back(index(x, 0, z));
            idx.push_back(index(x, NY_-1, z));
        }
    }
    for (int y = 0; y < NY_; ++y) {
        for (int x = 0; x < NX_; ++x) {
            idx.push_back(index(x, y, 0));
            idx.push_back(index(x, y, NZ_-1));
        }
    }
  /*
   for (int x = 0; x < NX_; ++x) {
            idx.push_back(index(x, 0, 0));
            idx.push_back(index(x, NY_-1, 0));
   }
    
   for (int y = 0; y < NY_; ++y) {
            idx.push_back(index(0, y, 0));
            idx.push_back(index(NX_-1, y, 0));
   }*/
   for (auto elem : idx) {
        std::cout << elem << " ";
    }
    std::cout <<"boundary" <<std::endl;
    for (int i = 0; i < idx.size(); ++i) {
        boundary_[idx[i]] = 1;
    }
    return boundary_;
}

int main () {
    Grid grid(3,3,3);
    std::cout << "total dimension: " << grid.dimension() << std::endl;
    std::cout << "spaceDime(): " << grid.spaceDim() << std::endl;
    std::cout << "boundary:\n";
    std::vector<int> bd;
    bd = grid.boundary();
    for (auto i : bd) {
        std::cout << i << "  ";
    }
    std::cout << std::endl;
}
