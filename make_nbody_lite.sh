#!/bin/bash    

# // Copyright Tom Donlon 2021 //
# //           RPI             //

#==========================================================
#This script is meant to be run whenever you want to make a new
#version of Milkyway@home nbody_lite. This script optionally (re)builds
#the MilkyWay@home Nbody binaries, and then (re)packages the folder
#mwah_nbody_lite, which will contain the files that should
#be released.
#==========================================================

#==========================================================
#$rebuild will build the binaries. Only needs to be done if
#milkywayathome has not been built yet, or if files have been
#changed or moved around since building.

#$package will create the mwah_nbody_lite folder (ready for release)
#prepared with all the binaries, documentation, and scripts that
#are needed to run MilkyWay@home Nbody.
#==========================================================

rebuild=true
package=true

#should always end up in the parent folder of the milkywayathome_client folder
cd ..

if $rebuild
then
    echo "Nbody Lite: Rebuilding MilkyWay@home Nbody binaries"
    rm -r build
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DNBODY_DEV_OPTIONS=ON -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DSEPARATION=OFF -DDOUBLEPREC=ON -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF -DAUTODIFF=OFF -DAUTODIFF_LOG=OFF ../milkywayathome_client/
    make -j 
    cd ..
fi

if $package
then
    #(re)build directory structure
    echo "Nbody Lite: Creating Nbody Lite directory structure"
    rm -r mwah_nbody_lite
    mkdir mwah_nbody_lite
    cd ./mwah_nbody_lite
    mkdir bin
    cd ..

    echo "Nbody Lite: Moving files to Nbody Lite"
    #get binary file
    cp ./build/bin/milkyway_nbody ./mwah_nbody_lite/bin

    #get .lua file and GUI .lua file
    cp ./milkywayathome_client/nbody/sample_workunits/settings.lua ./mwah_nbody_lite/bin
    cp ./milkywayathome_client/nbody/sample_workunits/settings_gui.lua ./mwah_nbody_lite/bin

    #get example manual bodies file
    cp ./milkywayathome_client/nbody/sample_workunits/manual_bodies_example.in ./mwah_nbody_lite/bin

    #get mwah_lite documentation, run.sh, and GUI files
    cp -r ./milkywayathome_client/lite/. ./mwah_nbody_lite
  
    echo "Nbody Lite: Nbody Lite successfully built"
fi
    
