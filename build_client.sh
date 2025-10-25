#!/bin/bash

cd "$(dirname "$0")"
PathToMilkyWayAtHomeClientDirectory="$(pwd)"
echo "Path to milkywayathome_client directory: $PathToMilkyWayAtHomeClientDirectory"

mkdir -p output
mkdir -p input

rm -r build
mkdir build
cd build
cmake \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS_DEBUG="-g -O0" \
  -DCMAKE_CXX_FLAGS_DEBUG="-g -O0" \
  -DNBODY_DEV_OPTIONS=ON \
  -DNBODY_GL=OFF \
  -DBOINC_APPLICATION=OFF \
  -DSEPARATION=OFF \
  -DDOUBLEPREC=ON \
  -DNBODY_OPENMP=ON \
  -DNBODY_OPENCL=OFF \
  "$PathToMilkyWayAtHomeClientDirectory/"
make -j
