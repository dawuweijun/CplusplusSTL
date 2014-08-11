# - Find OPM lattice library
#
# Defines the following variables:
#   opm-lattice_INCLUDE_DIRS    Directory of header files
#   opm-lattice_LIBRARIES       Directory of shared object files
#   opm-lattice_DEFINITIONS     Defines that must be set to compile
#   opm-lattice_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_OPM_LATTICE            Binary value to use in config.h

# Copyright (C) 2013 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

include (opm-lattice-prereqs)
include (OpmPackage)
find_opm_package (
  # module name
  "opm-lattice"

  # dependencies
  "${opm-lattice_DEPS}"

  # header to search for
  "opm/lattice/DataMap.hpp"

  # library to search for
  "opmlattice"

  # defines to be added to compilations
  ""

  # test program
"#include <opm/lattice/PolymerState.hpp>
int main (void) {
  DataMap s;
  return 0;
}
"
  # config variables
  "${opm-lattice_CONFIG_VAR}"
  )

#debug_find_vars ("opm-lattice")
