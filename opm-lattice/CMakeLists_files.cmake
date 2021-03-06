# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
#	MAIN_SOURCE_FILES     List of compilation units which will be included in
#	                      the library. If it isn't on this list, it won't be
#	                      part of the library. Please try to keep it sorted to
#	                      maintain sanity.
#
#	TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
#	TEST_DATA_FILES       Files from the source three that should be made
#	                      available in the corresponding location in the build
#	                      tree in order to run tests there.
#
#	EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#	                      build, but which is not part of the library nor is
#	                      run as tests.
#
#	PUBLIC_HEADER_FILES   List of public header files that should be
#	                      distributed together with the library. The source
#	                      files can of course include other files than these;
#	                      you should only add to this list if the *user* of
#	                      the library needs it.

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	opm/lattice/GridManager.cpp
	opm/lattice/LatticeBoltzmannModule.cpp
	opm/lattice/FluidProperties.cpp
	opm/lattice/LatticeBoltzmannSolver.cpp
	opm/lattice/LatticeBoltzmannSolverOutput.cpp
	opm/lattice/TwoPhaseLatticeBoltzmannSimulator.cpp
	opm/lattice/utility/writeVtkData.cpp
	opm/lattice/utility/StopWatch.cpp
	opm/lattice/utility/SimulatorState.cpp
	opm/lattice/utility/SimulatorTimer.cpp
	opm/lattice/utility/initState.cpp
	opm/lattice/utility/updateState.cpp
	)

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
	examples/test_grid.cpp
	examples/test_module.cpp
	examples/test_fluid.cpp
	)

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
	)

# originally generated with the command:
# find examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/sim_lbm.cpp
	)

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
	examples/sim_lbm.cpp
	#examples/hello.cpp
	)

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
	opm/lattice/GridManager.hpp
	opm/lattice/LatticeBoltzmannModule.hpp
	opm/lattice/FluidProperties.hpp
	opm/lattice/LatticeBoltzmannSolver.hpp
	opm/lattice/LatticeBoltzmannSolverOutput.hpp
	opm/lattice/TwoPhaseLatticeBoltzmannSimulator.hpp
	opm/lattice/utility/writeVtkData.hpp
	opm/lattice/utility/DataMap.hpp
	opm/lattice/utility/StopWatch.hpp
	opm/lattice/utility/SimulatorState.hpp
	opm/lattice/utility/SimulatorTimer.hpp
	opm/lattice/utility/initState.hpp
	opm/lattice/utility/updateState.hpp
	)
