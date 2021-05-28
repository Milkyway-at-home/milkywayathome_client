 // Copyright Tom Donlon 2021 //
 //           RPI             //

======================================================
MilkyWay@home Lite
======================================================

This application is meant to be a simpler approach to the MilkyWay@home
Nbody application. This is a way to install and use 
the Nbody project without having to go through and build the
entire Nbody project on your own, which has historically been a 
rather painful rite of passage for new users. 

Currently, MilkyWay@home Nbody Lite is only tested on and supported
for machines running the most recent version of Ubuntu OS (In theory 
it works on older versions of Ubuntu and different linux OS's, but this 
has not been tested). Windows and MacOS are not currently supported.  

======================================================
Installation
======================================================

 - Download the most recent nbody_lite release from 
https://github.com/Milkyway-at-home/milkywayathome_client/releases, 
or from milkyway.cs.rpi.edu (this second option may or may not be 
working as of the time you read this)

 - Extract the downloaded archive

======================================================
Usage
======================================================

 - Change settings in ./bin/settings.lua 
(optional, which settings you change depends on what type
of simulation you are running, and the parameters for that
specific simulation. The parameters that you can change
are explained in the settings.lua file.)

 - Open Terminal in the mwah_nbody_lite folder and run

./run.sh

This will start your simulation. It may be necessary to run

chmod +x run.sh

once before executing the run.sh file if you get an error that
you do not have the necessary permissions to execute run.sh. 

======================================================
What does Nbody Lite do?
======================================================

Nbody Lite has several components:

(i) - The software generates a dwarf galaxy with parameters that are 
given by the user. 

(ii) - The software integrates the dwarf galaxy particles backwards in
time in an external Milky Way potential.

(iii) - The software integrates the dwarf galaxy particles, and any
other manually-specified particles, forwards in time in an external 
Milky Way potential, including self gravity between the particles. 



The backwards and forward orbit integrations are performed so that a user
can initialize a simulation using either the start or end point of a dwarf
galaxy progenitor. If you want to specify the starting point of the dwarf, 
then set the reverse time to a very small number; if you want to specify 
the end point, then set the reverse time and forward time to be equal.

It is not necessary to run all of these steps: for example, 
one can choose not to run the dwarf galaxy particles backwards in time, one 
can avoid generating a dwarf galaxy at all in favor of using manual bodies, 
or the simulation can be run without an external potential. These settings can
all be changed in the settings.lua file. 

Nbody Lite now includes an optional LMC potential in simulations, 
which takes into account the gravity of the LMC on particles for reverse
and forward orbit integration. Nbody Lite also considers the reflex motion 
of the Milky Way from the LMC, as well as dynamical friction in the orbit
of the LMC. 

======================================================
Manual Body Input
======================================================

Nbody Lite supports the input of a list of bodies with manually
determined positions, velocities, and masses. These files must 
be in a specific format (tab separated value, or TSV) with a 
simple header. 

An example file with the correct formatting is provided in the
./bin folder as "manual_bodies_example.in". 

======================================================
FAQ
======================================================

Q: What is a structure mass?

A: A structure mass is a unit of mass used internally by 
MilkyWay@home that simplifies internal calculations. One
structure mass is equal to 22288.47 Solar masses.

