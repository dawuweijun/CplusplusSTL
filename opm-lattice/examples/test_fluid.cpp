#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GridTest

#include <boost/test/unit_test.hpp>
#include <opm/lattice/FluidProperties.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>

#include <iostream>

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(CreateGrid)
{
    LatticeBoltzmannModule module;
    FluidProperties red(module, 1.0, 1.0, 0.1667, 0.1);
    BOOST_CHECK_EQUAL( 1.0, red.rho());
    BOOST_CHECK_EQUAL( 1.0, red.tau());
    BOOST_CHECK_EQUAL( 0.1667, red.mu());
    BOOST_CHECK_EQUAL( 0.1, red.velmax());
}

BOOST_AUTO_TEST_SUITE_END()
