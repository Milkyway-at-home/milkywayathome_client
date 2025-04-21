#!/bin/bash

run=true
run_compare=false
compare_only=false
get_flag_list=false

cd "$(dirname "$0")"
PathToMilkyWayAtHomeClientDirectory="$(pwd)"
echo "Path to milkywayathome_client directory: $PathToMilkyWayAtHomeClientDirectory"

cd build/bin

if $run 
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientDirectory/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientDirectory/output/output.out \
    -z $PathToMilkyWayAtHomeClientDirectory/output/output.hist \
    -n 8 -b -w 1 -P -e 56551651 \
    -i 3.6333 1.0 0.1812 0.1828 1.223 0.0126 449865.888 \
    
fi

if $run_compare
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientDirectory/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientDirectory/output/output.out \
    -z $PathToMilkyWayAtHomeClientDirectory/output/output.hist \
    -h $PathToMilkyWayAtHomeClientDirectory/input/in3.hist \
    -n 8 -b -w 1 -P -e 44331451 \
    -p 3.6333 1.0 0.1812 0.1828 1.223 0.0126 \

fi

#SMU = 222,288.47 SOLAR MASSES

#OPTIONS:
#-s -> compare using only emd and cost component
#-S -> use emd, cost, beta dispersion
#-V -> use emd, cost, velocity dispersion
#-D -> use emd, cost, beta dispersion and velocity dispersion
if $compare_only 
then
    ./milkyway_nbody \
    -h $PathToMilkyWayAtHomeClientDirectory/input/in3.hist \
    -A $PathToMilkyWayAtHomeClientDirectory/output/output.hist \

fi

# if you run:

if $get_flag_list
then
    ./milkyway_nbody --help
# it will show you what all the flags mean
fi
