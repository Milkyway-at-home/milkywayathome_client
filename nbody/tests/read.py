import sys
args = sys.argv
f = open(args[2],'w')
sum_x = [0,0,0]
for line in open(args[1]):
		cats = line.split(',')
		if(len(cats)>8):
			x = float(cats[2])
			y = float(cats[3])
			z = float(cats[4])
			f.write(cats[2] + '\t' + cats[3] +'\t' +cats[4] +'\n')
			#ff.write(str(((x-sum_x[0])**2 + (y-sum_x[1])**2 + (z-sum_x[2])**2)**.5)+"\n")
		elif(len(cats) == 7):
			sum_x[0] = float(cats[0].split('=')[1])
			sum_x[1] = float(cats[1])
			sum_x[2] = float(cats[2])
print(sum_x)
