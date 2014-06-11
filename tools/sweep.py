#!/usr/bin/python

import collections
import subprocess
import sys
import os
import re

# Ugly script to run a parameter sweep for nbody and one of the existing workunits

# Print debug statements?
debug = False

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
luaFile = "../nbody/sample_workunits/EMD_10k_isotropic2.lua"
sweepName = sys.argv[1]
histDir = "./" + sweepName + "/hist"

# Create the directory if it doesn't exist
if not os.path.isdir(histDir):
    os.makedirs(histDir)

# Set up parameters to sweep
# Passed in the order they are added to the dictionary, this must 
# match the order in the lua file
params = collections.OrderedDict()

#       NAME                MIN   DEFAULT  MAX
params["forwardTime"]    = [1.9,  2,       2.1]
params["reverseRatio"]   = [0.95, 1,       1.05]
params["radius"]         = [0.9,  1,       1.1]
params["lightRadius"]    = [0.4,  0.5,     0.6]
params["mass"]           = [9,    10,      11]
params["lightMassRatio"] = [0.4,  0.5,     0.6]

# Return the default params for the run
def getDefaultParams(params):
    command = []
    for key in params:
        command += [ str(params[key][1]) ]
    return command

# Generate reference histogram
print "Generating reference histogram"
refHist = histDir + "/ref.hist"
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
    start = params[key][0]
    end = params[key][2]
    inc = (end - start) / float(steps)
    print start, ':', inc, ':', end

    f = open("./" + sweepName + "/" + key + '.txt', 'w')
    f.write("#" + str(getDefaultParams(params)) + "\n")

    # Sweep over the range of values
    for i in range(steps + 1):

        currentValue = start + i * inc
        f.write(str(currentValue) + "\t")

        # Save the output to the hist directory 
        histFile = histDir + "/" + key + "_" + str(currentValue) + ".hist"
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
        if( likelihood < -9999990 ):
            f.write("#")

        f.write(str(likelihood) + "\n")

    index += 1
    f.close()
