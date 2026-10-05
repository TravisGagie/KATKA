#!/bin/bash
# Builds Cliffy (pfp_doc / pfp_doc64) on this machine. Log: build_cliffy.log
cd "$(dirname "$0")/cliffy" || exit 1
export CMAKE_POLICY_VERSION_MINIMUM=3.5      # its dependencies ask for CMake versions older than 3.5
rm -rf build; mkdir -p build && cd build
# GCC 15 fixes go through a compiler wrapper, because some dependencies overwrite CMAKE_CXX_FLAGS
{ cmake -DCMAKE_CXX_COMPILER="$(cd ../.. && pwd)/gxx-compat" .. && make -j8 install; } > ../../build_cliffy.log 2>&1
st=$?
ls -la ../build/cliffy >> ../../build_cliffy.log 2>&1
echo "exit status $st" >> ../../build_cliffy.log
tail -5 ../../build_cliffy.log
