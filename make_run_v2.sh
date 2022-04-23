#!/bin/bash
#/* Copyright (c) 2016 Siddhartha Shelton */

rebuild=true
run=true
run_compare=false
compare_only=false
get_flag_list=false

PathToMilkyWayAtHomeClientFolder='INSERT PATHWAY HERE'

if $rebuild
then
    rm -r build
    mkdir build
    cd build
    cmake  -DCMAKE_BUILD_TYPE=Release -DNBODY_DEV_OPTIONS=ON -DNBODY_GL=OFF -DBOINC_APPLICATION=ON -DSEPARATION=OFF -DDOUBLEPREC=ON -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF -DAUTODIFF=OFF $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/
    make -j 
fi

cd /home/axiomomen/MilkyWay/build/bin

if $run 
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientFolder/test.out \
    -z $PathToMilkyWayAtHomeClientFolder/test.hist \
    -n 8 -b -w 1 -P -e 54231651\
    -i 4.0 1.0 0.2 0.2 12.0 0.2 \
    
fi

if $run_compare
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientFolder/output.out \
    -z $PathToMilkyWayAtHomeClientFolder/output.hist \
    -h $PathToMilkyWayAtHomeClientFolder/test.hist \
    -y $PathToMilkyWayAtHomeClientFolder/output.ad \
    -n 8 -b -w 1 -P -e 54231651\
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
    -h $PathToMilkyWayAtHomeClientFolder/test.hist \
    -S $PathToMilkyWayAtHomeClientFolder/output.hist \

fi


# if you run:

if $get_flag_list
then
    ./milkyway_nbody --help
# it will show you what all the flags mean
fi
