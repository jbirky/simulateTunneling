import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('classic')
rc('font', family='serif')
rc('figure', facecolor='w')
import os
import math
from math import sqrt, pi, exp, pow
import argparse

ALPHA = 1
BETA  = 2


def potential(x):

	pot = []
	for i in range(len(x)):

		# pot.append(1 / (ALPHA * sqrt(2*pi)) * ( exp(-pow((x[i] - BETA),2) / (2*pow(ALPHA,2))) + exp(-pow((x[i] + BETA),2) / (2*pow(ALPHA,2))) ))
		pot.append(ALPHA * pow(x[i], 4) - BETA * pow(x[i], 2) + pow(BETA, 2) / (4*ALPHA))

	return pot



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Specify plotting directory.')
	parser.add_argument("folder", action="store", type=str)
	args = parser.parse_args()


	# Plot expected values for each time step
	if 'expected' in args.folder:

		exp_files = os.listdir('expected/')

		for file in exp_files:
			vals = []
			with open('expected/' + file) as f:
				for line in f:
					vals.append(float(line))

			tstep = np.arange(0,len(vals),1)

			fname = file.split('.dat')[0]

			plt.figure(figsize=[12,6])
			plt.plot(tstep, vals)
			plt.xlabel('Timestep', fontsize=15)
			plt.ylabel(fname, fontsize=15)
			plt.savefig('plots/' + fname + '.png')
			# plt.show()
			plt.close()


		plt.figure(figsize=[12,6])
		for file in exp_files:
			vals = []
			with open('expected/' + file) as f:
				for line in f:
					vals.append(float(line))

			tstep = np.arange(0,len(vals),1)

			fname = file.split('.dat')[0]

			plt.plot(tstep, vals, label=fname, alpha=.8)
		plt.xlim(0,16)
		plt.xlabel('Timestep', fontsize=15)
		plt.ylabel('Expected Value', fontsize=15)
		plt.legend(loc='lower right')
		plt.title('Expected values vs. Time', fontsize=18)
		plt.savefig('plots/expected.png')
		# plt.show()
		plt.close()


	# Plot probability functions 
	if 'wave_prob' in args.folder:

		try:
			avg_eng_array = []
			with open('expected/avg_eng.dat') as f:
				for line in f:
						avg_eng_array.append(float(line))

			avg_eng_val = np.mean(avg_eng_array)
			show_eng = True
		except:
			show_eng = False

		prob_files = os.listdir('wave_prob/')

		for file in prob_files:
			vals = []
			with open('wave_prob/' + file) as f:
				for line in f:
					vals.append(float(line))
			vals = np.array(vals)

			xvals = np.linspace(-3, 3, len(vals))

			fname = file.split('.dat')[0]
			num = fname.split('phi_sq')[1]

			plt.figure(figsize=[12,8])
			plt.plot(xvals, potential(xvals), color='k', alpha=.8, label=r'$V(x) = \alpha x^4 - \beta x^2 + \frac{\beta^2}{4 \alpha}$')
			if show_eng == True:
				plt.plot(xvals, vals + avg_eng_val, label=r'$t_{s} = %s$'%(num), color='r')
				plt.axhline(y=avg_eng_val, linestyle='--', color='b', alpha=.4, label=r'$<E>$')
			else:
				plt.plot(xvals, vals, label=r'$t_{s} = %s$'%(num), color='r')
			plt.axvline(x=sqrt(BETA/(2*ALPHA)), linestyle='--', color='k', alpha=.4)
			plt.axvline(x=-sqrt(BETA/(2*ALPHA)), linestyle='--', color='k', alpha=.4)

			plt.ylim([-.2, 1.2])
			plt.xlabel(r'$x$', fontsize=15)
			plt.ylabel(r'$|\Psi(x,t)|^{2}$', fontsize=15)
			plt.legend(loc='upper right')
			plt.title('Wave Probability Amplitude', fontsize=18)
			plt.savefig('plots/phi/' + fname + '.png')
			plt.show()
			plt.close()

