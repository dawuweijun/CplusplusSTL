#!/bin/bash
# Download and build source code on RH5

	if [-d "build"]; then
		rm -rf build/
		else
		    mkdir build
		fi
    cd build
	cmake -DCMAKE_C_COMPILER="gcc48" -DCMAKE_CXX_COMPILER="gcc48" -DDUNE_ROOT="/project/res/x86_64_RH_5/share/opm/dune" -DSUPERLU_ROOT="/project/res/x86_64_RH_5/share/opm/SuperLU_4.0" -DBOOST_ROOT="/project/res/x86_64_RH_5/share/opm/boost" -DERT_ROOT="/project/res/x86_64_RH_5/share/ert/release/nightly/upstream/master" ../
