#!/bin/bash          

# // Copyright Tom Donlon 2021 //
# //           RPI             //

#==========================================================

#This script will run an Nbody simulation based on the given parameters
#below and in the for_developers.lua file (located in the ./bin folder)

#for -i, the parameters are:
# <evolveTime, revTime, rscale_l, light_r_ratio, mass_l, light_mass_ratio>
#optional parameters: 
# <... l, b, r, vx, vy, vz, manual_bodies.in> 
#the meaning of all of these parameters are provided in for_developers.lua
#NOTE: the manual_bodies file can be provided even the orbital parameters are not.

#valid examples:
# -i 1.0 1.0 0.3 0.2 45 0.1 
# -i 1.0 1.0 0.3 0.2 45 0.1 extra_particles.in
# -i 1.0 1.0 0.3 0.2 45 0.1 200.0 30.0 15.0 -50.0 100.0 -75.0 
# -i 1.0 1.0 0.3 0.2 45 0.1 200.0 30.0 15.0 -50.0 100.0 -75.0 extra_particles.in

#for info on other flags, you can run:
#cd ./bin
#./milkyway_nbody --help

#==========================================================

cd ./bin
./milkyway_nbody \
-f ./for_developers.lua \
-o output.out \
-n 8 -b -P \
-i 1.0 1.0 0.3 0.2 45 0.1 \

