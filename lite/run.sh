#!/bin/bash          

# // Copyright Tom Donlon 2021 //
# //           RPI             //

#==========================================================

#This script will run an Nbody simulation based on the given parameters
#in the settings.lua file (located in the ./bin folder)

#for info on other flags, you can run:
#cd ./bin
#./milkyway_nbody --help

#==========================================================

cd ./bin
./milkyway_nbody \
-f ./settings.lua \
-o output.out \
-n 8 -b -P \

