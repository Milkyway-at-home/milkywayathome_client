 // Copyright Tom Donlon 2021 //
 //           RPI             //

======================================================
MilkyWay@home Lite
======================================================

This application is meant to be a simpler approach to the MilkyWay@home
Nbody application. This is a way to install and use 
the Nbody project without having to go through and build the
entire Nbody project on your own, which has historically been a 
rather painful rite-of-passage for new users. 

Currently, MilkyWay@home Nbody Lite is only tested and supported
for machines running Ubuntu OS. (It works in theory on different 
linux OS's, but this has not been tested). Windows and MacOS are 
not currently supported. 

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
Manual Body Input
======================================================
Nbody Lite supports the input of a list of bodies with manually
determined positions, velocities, and masses. These files must 
be in a specific format (tab separated value, or TSV) with a 
simple header. 

An example file with the correct formatting is provided in the
./bin folder as "manual_bodies_example.in". 
 
