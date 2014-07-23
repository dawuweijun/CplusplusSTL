#ifndef GRIDMANAGER_HEADER_INCLUDED
#define GRIDMANAGER_HEADER_INCLUDED
/*
author: Ming Liu
email:  qilicun@outlook.com
*/
#include <vector>

class GridManager{
public:
    GridManager(int NX, int NY, int NZ);
    GridManager(int NX, int NY);
    ~GridManager();
    int dimension() const;
    int spaceDim() const;
    int index(const int x, const int y, const int z) const;
    const std::vector<int>& boundary() const { return boundary_; }
    std::vector<int>& boundary() {return boundary_; }
    const int NX() const { return NX_; }
    const int NY() const { return NY_; }
    const int NZ() const { return NZ_; }
     
private:
    int NX_;
    int NY_;
    int NZ_;
    std::vector<int> boundary_;
    void setBoundary();
};
#endif //GRID_HEADER_INCLUDED
