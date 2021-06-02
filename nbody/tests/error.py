import sys
import numpy as np
args = sys.argv
f1 = open(args[1])
f2 = open(args[2])
f3 = open(args[3],'w')
full = []
while(1):
	k = f1.readline()
	j = f2.readline()
	
	if not k:
		break
	d = np.sqrt(np.sum(np.square(np.abs(np.array(k.split(),dtype=np.float64) - np.array(j.split(),dtype=np.float64)))))
	full.append(d)
	f3.write(str(d))
	f3.write("\n")
full = np.array(full)
print(np.average(full), np.std(full))
	
