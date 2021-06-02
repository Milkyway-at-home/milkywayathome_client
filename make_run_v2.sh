#!/bin/bash          
#/* Copyright (c) 2016 Siddhartha Shelton */

rebuild=true
run=true
run_compare=false
compare_only=false
get_flag_list=false

if $rebuild
then
    cp /home/dylansheils/Desktop/Milkyway/build/bin/nbody_checkpoint /home/dylansheils/Desktop/temp
    rm -r build
    mkdir build
    cd build
    cmake  -DCMAKE_BUILD_TYPE=Debug -DNBODY_DEV_OPTIONS=ON -DEBUG=ON -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DSEPARATION=OFF -DDOUBLEPREC=ON -DNBODY_OPENMP=OFF -DNBODY_OPENCL=OFF /home/dylansheils/Desktop/Milkyway/milkywayathome_client
    make -j 
    cp /home/dylansheils/Desktop/temp/nbody_checkpoint /home/dylansheils/Desktop/Milkyway/build/bin/
    rm -r /home/dylansheils/Desktop/temp
    cd /home/dylansheils/Desktop/
    mkdir temp
fi

cp /home/dylansheils/Desktop/settings.lua /home/dylansheils/Desktop/Milkyway/build/bin
cd /home/dylansheils/Desktop/Milkyway/build/bin
if $run 
then

   ./milkyway_nbody -f ./settings.lua \
-o output.out \
-n 16 -w 1 -b -P \
    
fi


if $run_compare
then
    ./milkyway_nbody \
    -f /home/dylansheils/Desktop/Milkyway/EMD_v176_newdata2.lua \
    -o /home/dylansheils/Desktop/Milkyway/results/outputs/histogram_data/results_7854614814_datarun_6.out \
    -h /home/dylansheils/Desktop/Milkyway/results/inputs/data_hist_summer_2020_beta_disp3.hist \
    -z /home/dylansheils/Desktop/Milkyway/results/outputs/histogram_data/results_7854614814_datarun_6.hist \
    -n 8 -b -w 1 -P -e 34086709 \
    -i 3.64173 1.0 0.162541 0.219279 1.37384 0.0483034 \

fi

#SMU = 222,288.47 SOLAR MASSES

#-----------------------------------------------------------------------------------------------------------
#NEW SIMULATIONS
#-----------------------------------------------------------------------------------------------------------
#DATA 1 - 3.64106262096966, 1, 0.159881625801071, 0.218191882050774, 1.39560257998073, 0.05215907781124
#DATA 2 - 3.64069352302304, 1, 0.156612968990398, 0.211793921878332, 1.39262014782816, 0.0508805170139489
#DATA 3 - 3.64026514967517, 1, 0.161133327705519, 0.219702476263174, 1.3929082755872, 0.0534199987821891
#DATA 4 - 3.63974546780712, 1, 0.161674078463018, 0.219465492839004, 1.38409591322463, 0.0541437808727521
#DATA 5 - 3.63928853937945, 1, 0.162549532103122, 0.216249841039259, 1.39952502891379, 0.0482379795235213
#DATA 6 - 3.64172802417245, 1, 0.16254081320639, 0.219278553436109, 1.37384266510089, 0.04830341253167


#OPTIONS:
#-s -> compare using only emd and cost component
#-S -> use emd, cost, beta dispersion
#-V -> use emd, cost, velocity dispersion
#-D -> use emd, cost, beta dispersion and velocity dispersion
if $compare_only 
then
    ./milkyway_nbody_1.76_x86_64-pc-linux-gnu__mt \
    -h /home/dylansheils/Desktop/Milkyway/results/inputs/test_nbody100000.hist \
    -S /home/dylansheils/Desktop/Milkyway/results/outputs/resultstemp.hist \

fi





# if you run:

if $get_flag_list
then
    ./milkyway_nbody --help
# it will show you what all the flags mean
fi
