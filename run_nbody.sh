#!/bin/bash

run=false
run_compare=true
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
    -o $PathToMilkyWayAtHomeClientDirectory/output/output_equal_sidd.out \
    -z $PathToMilkyWayAtHomeClientDirectory/output/output_equal_sidd.hist \
    -n 11 -b -w 1 -P -e 54231651 \
    -i 4.0 1.0 0.2 0.2 12.0 0.2 \
    
fi

if $run_compare
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientDirectory/nbody/sample_workunits/for_developers_norm.lua \
    -o $PathToMilkyWayAtHomeClientDirectory/output/output_new.out \
    -z $PathToMilkyWayAtHomeClientDirectory/output/output_new.hist \
    -h $PathToMilkyWayAtHomeClientDirectory/output/output_equal_sidd.hist \
    -n 11 -b -w 1 -P -e 54231651 \
    -p 4.0 1.0 0.2 0.2 12.0 0.2 \

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
    -h $PathToMilkyWayAtHomeClientDirectory/inputs/test_nbody100000.hist \
    -S $PathToMilkyWayAtHomeClientDirectory/outputs/resultstemp.hist \

fi

# if you run:

if $get_flag_list
then
    ./milkyway_nbody --help
# it will show you what all the flags mean
fi
