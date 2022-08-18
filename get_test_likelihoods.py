########################################
# This code is used to recalculate the #
#   likelihoods for the model tests.   #
########################################


#Import necessary packages

import math
import os
import os.path as pth
import numpy as np

#Define lists

PathwayToClientParentDir = 'INSERT PATHWAY HERE'

model_list = ["model_1","model_2","model_3","model_4","model_5","model_5_bounds_test","model_6","model_7","model_8","model_9","model_ninkovic","model_triaxial","model_newhist1","model_newhist2","model_newhist3","model_LMC","model_bar","model_LMC_bar"]
body_list = ["100","1024","10000"]
seed_list = ["670828913", "886885833", "715144259", "430281807", "543966758"]

#model_list = ["model_LMC_bar"]
#body_list = ["100","1024"]
#seed_list = ["670828913", "886885833"]

empty_seed_array = []
empty_body_array = []
likelihoods = []
for s in range(len(seed_list)):
	empty_seed_array.append(' ')
for t in range(len(body_list)):
	exec("empty_body_array.append("+str(empty_seed_array)+")")
for u in range(len(model_list)):
	exec("likelihoods.append("+str(empty_body_array)+")")

lua_pth = PathwayToClientParentDir+'/milkywayathome_client/nbody/tests/RunTestUnits.lua'
tmp_pth = PathwayToClientParentDir+'/likelitemp.out'
result_pth = PathwayToClientParentDir+'/test_likelihoods.txt'

#Loop through all tests

for i in range(len(model_list)):
	for j in range(len(body_list)):
		for k in range(len(seed_list)):
			#Change seed in lua file
			print("REWRITING LUA FILE................................")
			temp = open('/home/axiomomen/MilkyWay/lua_temp.lua','w+')
			lua_file = open(lua_pth,'r')
			for line in lua_file:
				if (line.startswith('local testSeed =')):
					temp.write('local testSeed = testSeeds['+str(k+1)+']\n')
				else:
					temp.write(line)
			lua_file.close()
			temp.close()
			os.system('rm '+lua_pth)
			os.rename('/home/axiomomen/MilkyWay/lua_temp.lua',lua_pth)

			#Run test and get likelihood
			print("Running "+model_list[i]+"__"+body_list[j]+"_test: SEED = "+seed_list[k])
			check = False
			os.system("sh get_test_likelihoods.sh "+model_list[i]+" "+body_list[j]+" "+tmp_pth)
			tmp_file = open(tmp_pth, 'r')
			test_num = 0
			for line in tmp_file:
				if (line.startswith('test ')):
					test_num = int(line.replace('test ',''))
				prefix = str(test_num)+": "
				if (line.startswith(prefix+"<search_likelihood>")):
					string = line.replace(prefix+"<search_likelihood>-",'')
					like_score = float(string.replace('</search_likelihood>',''))
					print("Likelihood = "+str(like_score))
					likelihoods[i][j][k] = str(like_score)
					check = True
			if(not check):
				assert False    #NO LIKELIHOOD FOUND!
			tmp_file.close()
			os.system('rm '+tmp_pth)


#Write data
results = open(result_pth,'w+')
for i in range(len(model_list)):
	results.write('   ["'+model_list[i]+'"] = {\n')
	for j in range(len(body_list)):
		results.write('      ["'+body_list[j]+'"] = {\n')
		for k in range(len(seed_list)):
			if (k+1 == len(seed_list)):
				results.write('         ["'+seed_list[k]+'"] = '+str(likelihoods[i][j][k])+'\n')
			else:
				results.write('         ["'+seed_list[k]+'"] = '+str(likelihoods[i][j][k])+',\n')
		if (j+1 == len(body_list)):
			results.write('      }\n')
		else:
			results.write('      },\n\n')
	if (i+1 == len(model_list)):
		results.write('   }\n')
	else:
		results.write('   },\n\n')
results.close()			








