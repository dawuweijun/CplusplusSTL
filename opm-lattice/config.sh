#!/bin/bash
# Download and build source code on RH5

rm -rf build/

	if [-d "build"]; then
		rm -rf build/
		else
		    mkdir build
		fi
    cd build
#	cmake -DCMAKE_C_COMPILER="gcc48" -DCMAKE_CXX_COMPILER="g++48" -DBOOST_ROOT="/project/res/x86_64_RH_5/share/opm/boost" ../
	cmake -DCMAKE_C_COMPILER="gcc44" -DCMAKE_CXX_COMPILER="g++44" -DBOOST_ROOT="/project/eor/user/liuming/local" ../
