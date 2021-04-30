#!/usr/bin/python

import collections
import subprocess
import sys
import os
import re

# Ugly script to run a parameter sweep for nbody and one of the existing workunits

# Print debug statements?
debug = True

# Store all the data for the sweep the directory specificed by argv[1]
if(len(sys.argv) < 2):
    print("USAGE: sweep.py sweep_name")
    print("Data will be saved in ./sweep_name/ directory")
    print("Modify the script as you see fit")
    exit(1)

# Print warning
print("WARNING: This program may take several hours to run, consider using a utility such as screen if running over ssh")

# Define Constants
steps = 20
luaFile = "../nbody/sample_workunits/EMD_20k_isotropic2_custom_histogram.lua"
sweepName = sys.argv[1]
histDir = "./" + sweepName + "/hist/"

# Create the directory if it doesn't exist
if not os.path.isdir(histDir):
    os.makedirs(histDir)

# Set up parameters to sweep
# Passed in the order they are added to the dictionary, this must 
# match the order in the lua file
params = collections.OrderedDict()

#       NAME               STARTING VALUE
params["forwardTime"]    = 2
params["reverseRatio"]   = 1
params["radius"]         = 1
params["lightRadius"]    = 0.5
params["mass"]           = 10
params["lightMassRatio"] = 0.5

# Return the default params for the run
def getDefaultParams(params):
    command = []
    for key in params:
        command += [ str(params[key]) ]
    return command

# Generate reference histogram
print "Generating reference histogram"
refHist = histDir + "ref.hist"
command = ["milkyway_nbody","-e", "0", "-f", luaFile, "-z", refHist, "--disable-opencl", "-i"]
command += getDefaultParams(params)
if(debug):
    print " ".join(command)
subprocess.call(command)

# Loop over all the parameters we are sweeping
index = 0
for key in params:

    # Print info about sweep, calculate increment
    print 'Sweeping parameter:', key
    start = params[key] * float(0.9)
    end = params[key] * float(1.1)
    inc = (end - start) / float(steps)
    print start, ':', inc, ':', end

    f = open("./" + sweepName + "/" + key + '.txt', 'w')
    f.write("#" + str(getDefaultParams(params)) + "\n")

    # Sweep over the range of values
    for i in range(steps + 1):

        currentValue = start + i * inc
        f.write(str(currentValue) + "\t")

        # Save the output to the hist directory 
        histFile = histDir + key + "_" + str(currentValue) + ".hist"
        command = ["milkyway_nbody","-e", "0", "-f", luaFile, "-z", histFile, "--disable-opencl", "-i", "-h", refHist]

        # Generate the params to pass to the current run
        currentParams = getDefaultParams(params)
        currentParams[index] = str(currentValue)
        command += currentParams

        if(debug):
            print " ".join(command)

        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        
        if(debug):
            print "STDOUT:", stdout
            print "STDERR:", stderr
        
        match = re.search("<search_likelihood>(.*)</search_likelihood>", stderr)
        likelihood = match.group(1)
        print(str(currentValue) + ":" +  str(likelihood))
        

        # Do not plot worst case likelihoods
        if( float(likelihood) < -100000 ):
            f.write("#")

        f.write(str(likelihood) + "\n")

    index += 1
    f.close()
