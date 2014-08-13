#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GridTest

#include <boost/test/unit_test.hpp>

//#include <boost/test/floating_point_comparsion.hpp>

#include <opm/lattice/GridManager.hpp>

#include <iostream>

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(CreateGrid)
{
    const int NX = 4;
    const int NY = 4;
    const int NZ = 4;
    GridManager grid(NX, NY, NZ);
    BOOST_CHECK_EQUAL( NX, grid.NX());
    BOOST_CHECK_EQUAL( NY, grid.NY());
    BOOST_CHECK_EQUAL( NZ, grid.NZ());
    BOOST_CHECK_EQUAL( 3, grid.spaceDim());
    BOOST_CHECK_EQUAL( NX*NY*NZ, grid.dimension());
}

BOOST_AUTO_TEST_SUITE_END()
