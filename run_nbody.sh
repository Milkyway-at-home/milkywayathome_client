#!/bin/bash

# CHANGE THESE TO RUN A BASIC SIMULATION
# - Lua configuration file, Input files, Output files (lines 26-28, e.g.)

# RUN WITH DEBUG OPTIONS (stack trace with line numbers)
# - add "gdb --args " before your ./milkyway_nbody call (below)
# - when (gdb) pops up in terminal, type "run", then "y" when prompted
# - to trace your error, type "bt" after the program fails

run=true
run_compare=false
compare_only=false
get_flag_list=false

cd "$(dirname "$0")"
PathToMilkyWayAtHomeClientDirectory="$(pwd)"
echo "Path to milkywayathome_client directory: $PathToMilkyWayAtHomeClientDirectory"

timestamp=$(date +%Y%m%d_%H%M%S)
cd build/bin

#Note that the server rounds all parameters to 6 significant figures before sending them to the client
#So if you want to test the exact same parameters as the server, round to 6
if $run 
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientDirectory/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientDirectory/output/output_${timestamp}.out \
    -z $PathToMilkyWayAtHomeClientDirectory/output/output_${timestamp}.hist \
    -n 8 -w 1 -P -e 54231651 \
    -i 4.0 1.0 0.2 0.2 12.0 0.2 \
    
fi

if $run_compare
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientDirectory/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientDirectory/output/output_${timestamp}.out \
    -z $PathToMilkyWayAtHomeClientDirectory/output/output_${timestamp}.hist \
    -h $PathToMilkyWayAtHomeClientDirectory/input/input.hist \
    -n 8 -w 1 -P -e 54231651 \
    -p 4.0 1.0 0.2 0.2 12.0 0.2 \

fi

#SMU = 222,288.47 SOLAR MASSES

#OPTIONS:
#-s -> The histogram to compare. Will always compare using emd and cost component
#ADD ADDITIONAL FLAGS TO INCLUDE OTHER COMPARISONS
#-S -> Include beta dispersion in the comparison
#-V -> Include velocity dispersion in the comparison
#-B -> Include beta average in the comparison
#-Q -> Include line of sight velocity average in the comparison
#-D -> Include distance average in the comparison
#-U -> Include proper motion average in the comparison
#-L -> Include momentum in the comparison
#Note that average beta bins with be read in with the histogram given with -h, 
#while EMDRange will prioritze the range given in the histogram given with -s
if $compare_only 
then
    ./milkyway_nbody \
    -h $PathToMilkyWayAtHomeClientDirectory/input/input.hist \
    -s $PathToMilkyWayAtHomeClientDirectory/output/output.hist \

fi

# if you run:

if $get_flag_list
then
    ./milkyway_nbody --help
# it will show you what all the flags mean
fi
