#!/bin/bash
#/* Copyright (c) 2016 Siddhartha Shelton */

rebuild=false
run=false
run_compare=true
compare_only=false
get_flag_list=false

PathToMilkyWayAtHomeClientFolder='INSERT PATHWAY HERE'

if $rebuild
then
    rm -r build
    mkdir build
    cd build
    cmake  -DCMAKE_BUILD_TYPE=Release -DNBODY_DEV_OPTIONS=ON -DNBODY_GL=OFF -DBOINC_APPLICATION=ON -DSEPARATION=OFF -DDOUBLEPREC=ON -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/
    make -j 
fi

cd $PathToMilkyWayAtHomeClientFolder/build/bin

if $run 
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientFolder/results/inputs/test.out \
    -z $PathToMilkyWayAtHomeClientFolder/results/inputs/test.hist \
    -n 8 -b -w 1 -P -e 54231651 \
    -i 4.0 1.0 0.2 0.2 12.0 0.2 \
    
fi

if $run_compare
then
    ./milkyway_nbody \
    -f $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/sample_workunits/for_developers.lua \
    -o $PathToMilkyWayAtHomeClientFolder/results/outputs/histogram_data/output.out \
    -z $PathToMilkyWayAtHomeClientFolder/results/outputs/histogram_data/output.hist \
    -h $PathToMilkyWayAtHomeClientFolder/results/inputs/test.hist \
    -n 8 -b -w 1 -P -e 54231651 \
    -p 4.0 1.0 0.2 0.2 12.0 0.2 \

fi

#SMU = 222,288.47 SOLAR MASSES

#-----------------------------------------------------------------------------------------------------------
#1.76 OPTIMIZATIONS
#-----------------------------------------------------------------------------------------------------------
#DATA 1 - 3.64106262096966, 1, 0.159881625801071, 0.218191882050774, 1.39560257998073, 0.05215907781124
#DATA 2 - 3.64069352302304, 1, 0.156612968990398, 0.211793921878332, 1.39262014782816, 0.0508805170139489
#DATA 3 - 3.64026514967517, 1, 0.161133327705519, 0.219702476263174, 1.3929082755872, 0.0534199987821891
#DATA 4 - 3.63974546780712, 1, 0.161674078463018, 0.219465492839004, 1.38409591322463, 0.0541437808727521
#DATA 5 - 3.63928853937945, 1, 0.162549532103122, 0.216249841039259, 1.39952502891379, 0.0482379795235213
#DATA 6 - 3.64172802417245, 1, 0.16254081320639, 0.219278553436109, 1.37384266510089, 0.04830341253167
#DATA 7 - 3.64260394155576, 1, 0.160612762606858, 0.226057114089681, 1.35028221950791, 0.0587807840203055
#DATA 8 - 3.64131941608212, 1, 0.198202679607823, 0.295589684251824, 1.38172044956831, 0.0685817884919194
#DATA 9 - 3.63670928512826, 1, 0.192445000718284, 0.277168017788357, 1.37774704258012, 0.0650096344878362

#DATA 1 - 3.63535222803186, 1, 0.255949559307738, 0.345500150804994, 1.01666060326807, 0.04131780899306
#DATA 2 - 3.63445771619562, 1, 0.232665840250609, 0.254335348225735, 1.14616007319711, 0.0179472956620577     ***
#DATA 3 - 3.63329947776065, 1, 0.18121623259995, 0.182798519641869, 1.22250889107504, 0.0126170901673599      ***
#DATA 4 - 5.39211743041288, 1, 0.205989430632811, 0.131469783691291, 1.63492542526056, 0.144474632927037
#DATA 5 - 3.63337849690933, 1, 0.184161939447946, 0.1820116923855, 1.25090662174001, 0.0119267105962461       ***
#DATA 6 - 5.39563974247424, 1, 0.19871047747507, 0.13053729109133, 1.68176202390101, 0.150690744489736
#DATA 7 - 2.72671181827783, 1, 0.0968623606530016, 0.169724369461875, 0.55257187727682, 0.0112672850567859
#DATA 8 - 5.39624868562447, 1, 0.200471628505676, 0.129628510794245, 1.72840723764473, 0.155166316711468
#DATA 9 - 3.63167873245535, 1, 0.152934244619535, 0.16413882014093, 1.18472334100967, 0.0153673674324755      ***
#DATA 10- 3.64423965699323, 1, 0.225262235814408, 0.365341742131286, 0.965902951499858, 0.0667439709840524
#DATA 11- 3.63521737931889, 1, 0.252275846008539, 0.343163951842872, 1.03082343117794, 0.0425020532412654
#DATA 12- 3.63224533886196, 1, 0.219497624725169, 0.246650735071382, 1.14458244557394, 0.0199777503652921     ***
#DATA 13- 3.63529957226671, 1, 0.25131609694193, 0.34309180304589, 1.02389680127682, 0.0424752671679215

#RUN---l-------------b------------ra------------dec
#1   267.87414641   41.94005915   166.76762953  -13.71819358   
#2   264.0278193    44.39400418   165.73958557  -10.13004331
#3   270.38493762   40.53692443   167.70186807  -15.89005682
#4   //
#5   260.36064825   46.02149916   164.59175198  -7.29559970 
#6   //

#OPTIONS:
#-s -> compare using only emd and cost component
#-S -> use emd, cost, beta dispersion
#-V -> use emd, cost, velocity dispersion
#-D -> use emd, cost, beta dispersion and velocity dispersion
if $compare_only 
then
    ./milkyway_nbody \
    -h $PathToMilkyWayAtHomeClientFolder/results/inputs/test_nbody100000.hist \
    -S $PathToMilkyWayAtHomeClientFolder/results/outputs/resultstemp.hist \

fi

# if you run:

if $get_flag_list
then
    ./milkyway_nbody --help
# it will show you what all the flags mean
fi
