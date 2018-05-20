import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('classic')
rc('font', family='serif')
rc('figure', facecolor='w')
import os
import math
from math import sqrt, pi, exp
import argparse
import shutil

V_ALPHA = .02
BETA  = 2
XMIN = 2.5


def potential(x):

	pot = []
	for i in range(len(x)):

		# pot.append(ALPHA * x[i]**4 - BETA * x[i]**2 + BETA**2 / (4*ALPHA))
		pot.append(V_ALPHA * (x[i]**2 - XMIN**2)**2)

	return pot



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Specify plotting directory.')
	parser.add_argument("plot", action="store", type=str)
	# parser.add_argument("show", action="store", type=bool)
	args = parser.parse_args()


	# Plot expected values for each time step
	if 'expected' in args.plot:

		print('Plotting expected values...')

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

		plt.axhline(y=XMIN, linestyle='--', color='k', alpha=.4, label=r'$\pm x_{min}$')
		plt.axhline(y=-XMIN, linestyle='--', color='k', alpha=.4)

		plt.xlim(0, len(vals)-1)
		plt.ylim(-XMIN-.1, XMIN+.1)
		plt.xlabel('Timestep', fontsize=15)
		plt.ylabel('Expected Value', fontsize=15)
		plt.legend(loc='lower left', fontsize=12)
		plt.title('Expected values vs. Time', fontsize=18)
		plt.savefig('plots/expected.png')
		# plt.show()
		plt.close()


	# Plot probability functions 
	if 'wave_prob' in args.plot:

		print('Plotting wave amplitudes...')

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

			xvals = np.linspace(-6, 6, len(vals))

			fname = file.split('.dat')[0]
			num = fname.split('phi_sq')[1]

			plt.figure(figsize=[12,8])
			plt.plot(xvals, potential(xvals), color='k', alpha=.8, label=r'$V(x) = \alpha (x^2 - x_{min}^2)^2$')
			if show_eng == True:
				plt.plot(xvals, vals + avg_eng_val, label=r'$|\Psi(x,t_{s} = %s)|^{2}$'%(num), color='r')
				plt.axhline(y=avg_eng_val, linestyle='--', color='b', alpha=.4, label=r'$<E>=%s$'%(str(round(avg_eng_val,3))))
			else:
				plt.plot(xvals, vals, label=r'$t_{s} = %s$'%(num), color='r')
			plt.axvline(x=XMIN, linestyle='--', color='k', alpha=.4)
			plt.axvline(x=-XMIN, linestyle='--', color='k', alpha=.4)
			plt.text(-5.8, 1.13, r'$\alpha=%s$, $x_{min}=%s$'%(str(V_ALPHA), str(XMIN)))

			plt.ylim([-.1, 1.2])
			plt.xlabel(r'$x$', fontsize=15)
			plt.ylabel(r'$|\Psi(x,t)|^{2}$', fontsize=15)
			plt.legend(loc='lower right')
			plt.title('Wave Probability Amplitude', fontsize=18)
			plt.savefig('plots/phi/' + fname + '.png')
			# plt.show()
			plt.close()


	if 'all' in args.plot:

		exp_files = os.listdir('expected/')

		expected = []
		for file in exp_files:
			vals = []
			with open('expected/' + file) as f:
				for line in f:
					vals.append(float(line))
			expected.append(vals)


	if 'period' in args.plot:

		show = [0, 32, 64, 96, 128]
		# show = [0, 16, 32, 48, 64]

		print('Plotting %s samples in the period...'%(str(len(show))))

		try:
			avg_eng_array = []
			with open('expected/avg_eng.dat') as f:
				for line in f:
						avg_eng_array.append(float(line))

			avg_eng_val = np.mean(avg_eng_array)
			show_eng = True
		except:
			show_eng = False

		plt.figure(figsize=[12,8])

		for num in show:
			vals = []
			with open('wave_prob/phi_sq%s.dat'%(str(num))) as f:
				for line in f:
					vals.append(float(line))
			vals = np.array(vals)
			xvals = np.linspace(-6, 6, len(vals))
			
			plt.plot(xvals, vals + avg_eng_val, alpha=.8)
		
		plt.plot(xvals, potential(xvals), color='k', alpha=.8, label=r'$V(x) = \alpha (x^2 - x_{min}^2)^2$')
		if show_eng == True:		
			plt.axhline(y=avg_eng_val, linestyle='--', color='b', alpha=.4, label=r'$<E>=%s$'%(str(round(avg_eng_val,3))))
		plt.axvline(x=XMIN, linestyle='--', color='k', alpha=.4)
		plt.axvline(x=-XMIN, linestyle='--', color='k', alpha=.4)
		plt.text(-5.8, 1.13, r'$\alpha=%s$, $x_{min}=%s$'%(str(V_ALPHA), str(XMIN)))

		plt.ylim([-.1, 1.2])
		plt.xlabel(r'$x$', fontsize=15)
		plt.ylabel(r'$|\Psi(x,t)|^{2}$', fontsize=15)
		plt.legend(loc='lower right')
		plt.title('Wave Probability Amplitude', fontsize=18)
		plt.savefig('plots/period_sample.png')
		# plt.show()
		plt.close()


	if 'fourier' in args.plot:

		from scipy.fftpack import fft

		array = []
		with open('expected/avg_pos.dat') as f:
			for line in f:
				array.append(float(line))

		ft = fft(array)
		tstep = np.arange(0, len(ft), 1)
		# tstep = np.linspace(0, 2*math.pi/128, 64)

		xmax = np.where(np.array(ft) == max(ft))[0][0]
		max_val = ft[xmax]

		plt.figure(figsize=[12,8])
		plt.plot(tstep, ft)
		plt.axvline(x=xmax, label=r'$\omega_{max}=%s$'%(str(max_val)), color='r', alpha=.7)

		plt.xlim(0, 30)
		plt.legend(loc='upper right')
		plt.savefig('plots/pos_fourier.png')
		# plt.show()
		plt.close()


	if 'animation' in args.plot:

		l = len(os.listdir("plots/phi/")) - 1

		print("Converting {} images to gif...".format(l))

		os.system("convert -delay 5 $(for i in $(seq 0 1 %s); do echo plots/phi/phi_sq${i}.png; done) \
			-loop 0  plots/wave_evolution.gif"%(str(l)))