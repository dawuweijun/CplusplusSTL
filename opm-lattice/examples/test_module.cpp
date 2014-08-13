#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ModuleTest

#include <boost/test/unit_test.hpp>
#include <opm/lattice/LatticeBoltzmannModule.hpp>

#include <iostream>

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(CreateModule)
{
    LatticeBoltzmannModule module;
    std::cout << "Module name: " << module.name() << std::endl;
    BOOST_CHECK_EQUAL( 19, module.numDirection());
}

BOOST_AUTO_TEST_SUITE_END()
