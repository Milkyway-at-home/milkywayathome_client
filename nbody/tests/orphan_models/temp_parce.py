#Import necessary packages
import matplotlib.pyplot as plt

#Input hist file

g = open('orphan_model_histogram', 'r')

#Parce data from histogram file into arrays

binnum = 0
matterhist = []
lambd = []
counts = []
error = []

for line in g:
	hs = line.split(' ')

	typ = int(hs[0])
	lamval = float(hs[1])
	countval = float(hs[2])
	errorval = float(hs[3])
		
	matterhist.append(typ)
	lambd.append(lamval)
	counts.append(countval)
	error.append(errorval)

	binnum = binnum + 1
	skipline = True
	#print('Lambda Bin ', binnum, ' logged.')

print('ALL LAMBDA BINS LOGGED.')

#Close imput file

g.close()

#Generate Lambda-Beta Histograms

plt.figure(1)

plt.bar(lambd, counts, width=3)		#1D LAMBDA PLOT
plt.ylabel('Normalized Counts')
plt.xlabel(r'$\Lambda$' + ' (degrees)')
plt.axis([50,-50,0,0.35])
plt.title('Distribution of Orphan Stream')

#Output histograms

plt.show()
