# Emily Crook 2020
# This code tests the MainStruct / NBodyHistogram set up
# The six histograms are: normalized counts, beta disp, vel disp,
# beta average, velocity average, and distance
# It also tests the simultaneous fitting of the orbit parameters
# and backwards compatibility

# IMPORTANT: CHANGE THE PATHWAYS IN THE SHELL SCRIPTS

import os
import os.path as pth

# possible error messages (not a catch-all but these are common ones)
error1 = 'Likelihood was NAN. Returning worst case. \n'
error2 = 'Segmentation fault (core dumped)\n'

# default values for being compared
params1 = "4.0 1.0 0.2 0.2 12 0.2 52.5 28.6 -156 79 107"
params2 = "4.0 1.0 0.2 0.2 12 0.2 52.5 28.6 0 79 107"
params3 = "4.0 1.0 0.2 0.2 12 0.2" # original input parameters

# run two tests with different parameters
# (both using the new histogram/likelihood)
# and make sure the likelihood goes down
os.system('sh newhiststruct_run.sh '+params1)
os.system('sh newhiststruct_compare.sh '+params2)

if(pth.isfile("likelihood.out")):
	f = open('/likelihood.out', 'r')
	string = f.readline()
	string1 = string.replace("<search_likelihood>", "")
	string2 = string1.replace("</search_likelihood>", "")

	if(string2 == error1 or string2 == error2):
		print("Error message received, first test failed. Error message: ", string2, "\n")
	else:
		likelihood = float(string2)
		if(likelihood < 0):
			print("First test successful, likelihood: ", likelihood)
		else:
			print("First test failed, identical histograms produced")
else:
	print("likelihood file not produced")

# delete generated files?