#!/bin/bash

cd "$(dirname "$0")"
PathToMilkyWayAtHomeClientDirectory="$(pwd)"
echo "Path to milkywayathome_client directory: $PathToMilkyWayAtHomeClientDirectory"

rm -r build
mkdir build
cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DNBODY_DEV_OPTIONS=ON -DNBODY_GL=OFF -DBOINC_APPLICATION=ON -DSEPARATION=OFF -DDOUBLEPREC=ON -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF $PathToMilkyWayAtHomeClientDirectory/
make -d
