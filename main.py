from numpy import array, linspace, log, exp, seterr
from math import log2
from numpy.random import power, lognormal, random, randint
import matplotlib.pyplot as plt
from matplotlib import style
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

style.use('classic')
seterr('ignore')
		
def imf(model, n):
	'''Initial mass function'''
	i = 0
	s = []
	M0 = 1.0
	M = [1, 10]
	beta = 2.35
	if model == 'A':
		#Uniform star mass
		return [M0 for _ in range(n)]
	elif model == 'B':
		#Salpeter mass distribution
		while i!=n:
			rnd = random()
			sM=(((1.0-rnd)/(M[0]**(beta-1.0))) + (rnd/(M[1]**(beta-1.0))))**((-1.0)/(beta-1.0))
			s.append(sM)
			if 10>sM>1:
				i+=1
		return s
	elif model == 'C':
		#Lognormal mass distribution
		for _ in range(n):
			sM = lognormal()+1
			s.append(sM)
			if 10>sM>1:
				i+=1
		return s 
		
def merge(stars):
	i1 = randint(0, len(stars))
	i2 = randint(0, len(stars))
	if i1!=i2:
		stars[i1] += stars[i2] 
		del stars[i2]
	else:
		merge(stars)


def main():
	print('-------------------starMerging-------------------')
	print('----------------------v0.0.1---------------------')
	print('You can choose initial mass function:')
	print('A - uniform IMF\nB - Salpeter IMF\nC - lognormal IMF')
	print('Press "Enter" to choose default parameters:')
	model = input('Choose mass function (default "B"):') or 'B'
	n = input('Enter number of stars (default 100000):') or 100000; n = int(n);
	p = input('Enter percent of merging (default 0.95):') or 0.95; p = float(p);
	print("Calculating... it may take few minutes, it's ok :)")
	L = int(p*n)

	stars = imf(model, n)

	def func(x, a, b, c):
		if model in '':
			return a*x**(-b)+c
		elif model in 'ABC':
			return a*exp(-b*x)+c

	if model == 'A':
		nbins = 10
	else:
		nbins = int(log2(n)+1)
	
	bounds = [1, 20]

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	plt.hist(stars, bins = [i for i in linspace(bounds[0], bounds[1], nbins)], label = 'IMF', color = '#0092CC', alpha = 0.7)
	plt.title('Initial and result mass distributions ({} stars)'.format(n))
	plt.xlabel(r'Mass, $M_\odot$ ')
	plt.ylabel(r'Number of stars, N')
	plt.legend()
	plt.savefig('imf.png')
	
	#Stars stars merging
	for _ in range(L):
		merge(stars)
	plt.hist(stars, bins = [i for i in linspace(bounds[0], bounds[1], nbins)], label = 'RMF', color = '#FF3333', alpha = 0.7)
	plt.legend()
	plt.savefig('compared.png')
	plt.clf()
	bounds = [int(min(stars)), int(max(stars)) if max(stars)<150 else 150]
	vals, bins, patches = plt.hist(stars, bins = [i for i in linspace(bounds[0], bounds[1], nbins)], color = '#0092CC')
	newBins = array([(bins[i]+bins[i+1])/2 for i in range(0,len(bins)-1)])
	
	try:
		popt, pcov = curve_fit(func, newBins, vals, maxfev = L)
		funcBins = linspace(bounds[0],bounds[1],100)
		plt.plot(funcBins, func(funcBins, *popt),
					label='Approximation: $ae^{-%1.3fM}+c$'%popt[1], color = '#FF3333')
		plt.legend()
	except:
		print("Approximation not found :(")

	plt.plot(newBins, vals, 'o', color = '#FF3333')
	plt.title('Mass distribution after {} ({} \%) merging'.format(L,100*p))
	plt.ylim(0,max(vals)+0.2*max(vals))
	plt.xlabel(r'Mass, $M_\odot$ ')
	plt.ylabel(r'Number of stars, N')
	
	plt.savefig('rmf.png')
	plt.show()
	print('Finish! If all was ok - you can find rmf.png and imf.png in work directory.')

main()

#imf('C', n)
#y = imf('B', n)
#print(choice(y))
#print(len(y))
#plt.ylim(0,1200)
#plt.hist(y, bins = [i for i in linspace(1, 10,1000)])
#plt.savefig('fig2.eps')
#y = imf('C', n)
#print(len(y))
#plt.ylim(0,1200)
#plt.hist(y, bins = [i for i in linspace(1, 10,1000)])
#plt.savefig('fig1.eps')
#plt.show()

#plt.show()

