import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from itertools import product
from scipy.interpolate import interp2d
import numpy as np
import os

def onebody_plot2 ():

	dir = os.path.abspath('..' + '/plots/')

	fig = plt.figure(frameon=False)

	ax = fig.add_subplot(111)

	ax.set_xlabel('${N}_{max}$', fontsize = 22)
	ax.set_ylabel('Energy [MeV]', fontsize = 22)

	fig.suptitle('2 particles - trap $\hbar \omega$ = 10 MeV', fontsize = 22)

	plot_name = 'onebody_results2'

	methods  = ['k-', 'r-',  'g-', 'b-', 'y-', 'm-', 'c-', 'r--', 'g--', 'b--']

	nmax = [2, 4, 6, 8, 10]

	hw10 = [30.00, 30.00, 30.00, 30.00, 30.00]
	hw15 = [30.22, 30.02, 30.00, 30.00, 30.00]
	hw20 = [31.48, 30.28, 30.05, 30.01, 30.00]
	hw25 = [33.74, 31.08, 30.30, 30.08, 30.02]
	hw30 = [36.67, 32.42, 30.87, 30.30, 30.10]
	hw35 = [40.03, 34.20, 31.79, 30.75, 30.30]
	hw40 = [43.68, 36.30, 33.00, 31.42, 30.66]
	hw45 = [47.52, 38.64, 34.45, 32.31, 31.19]
	hw50 = [51.51, 41.15, 36.09, 33.38, 31.88]
	hw60 = [59.79, 46.57, 39.80, 35.95, 33.65]

	plt.subplot(111).plot(nmax, hw60, methods[9], label = '$\hbar \omega = 60$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw50, methods[8], label = '$\hbar \omega = 50$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw45, methods[7], label = '$\hbar \omega = 45$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw40, methods[6], label = '$\hbar \omega = 40$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw35, methods[5], label = '$\hbar \omega = 35$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw30, methods[4], label = '$\hbar \omega = 30$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw25, methods[3], label = '$\hbar \omega = 25$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw20, methods[2], label = '$\hbar \omega = 20$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw15, methods[1], label = '$\hbar \omega = 15$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw10, methods[0], label = '$\hbar \omega = 10$ MeV', linewidth=2)

	ax.legend(loc = 'best', frameon=False, fontsize = 18)

	fig.set_size_inches(8,8)

	ax.xaxis.set_tick_params(labelsize=20)
	ax.yaxis.set_tick_params(labelsize=20)

	plt.savefig('plots/' + plot_name + '.pdf', format='pdf')


def onebody_plot4 ():

	dir = os.path.abspath('..' + '/plots/')

	fig = plt.figure(frameon=False)

	ax = fig.add_subplot(111)

	ax.set_xlabel('${N}_{max}$', fontsize = 22)
	ax.set_ylabel('Energy [MeV]', fontsize = 22)

	fig.suptitle('4 particles - trap $\hbar \omega$ = 10 MeV', fontsize = 22)

	plot_name = 'onebody_results4'

	methods  = ['k-', 'r-',  'g-', 'b-', 'y-', 'm-', 'c-', 'r--', 'g--', 'b--']

	nmax = [2, 4, 6, 8, 10]

	hw10 = [80.00, 80.00, 80.00, 80.00, 80.00]
	hw15 = [84.38, 80.49, 80.00, 80.00, 80.00]
	hw20 = [93.98, 83.32, 80.76, 80.15, 80.03]
	hw25 = [106.2, 88.49, 82.84, 80.91, 80.27]
	hw30 = [120.0, 95.34, 86.30, 82.57, 81.01]
	hw35 = [134.7, 103.3, 90.86, 85.12, 82.38]
	hw40 = [149.9, 112.1, 96.25, 88.46, 84.40]
	hw45 = [165.6, 121.4, 102.3, 92.42, 86.99]
	hw50 = [181.5, 131.2, 108.8, 96.89, 90.06]
	hw60 = [214.0, 151.5, 122.8, 106.9, 97.35]

	plt.subplot(111).plot(nmax, hw60, methods[9], label = '$\hbar \omega = 60$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw50, methods[8], label = '$\hbar \omega = 50$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw45, methods[7], label = '$\hbar \omega = 45$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw40, methods[6], label = '$\hbar \omega = 40$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw35, methods[5], label = '$\hbar \omega = 35$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw30, methods[4], label = '$\hbar \omega = 30$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw25, methods[3], label = '$\hbar \omega = 25$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw20, methods[2], label = '$\hbar \omega = 20$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw15, methods[1], label = '$\hbar \omega = 15$ MeV', linewidth=2)
	plt.subplot(111).plot(nmax, hw10, methods[0], label = '$\hbar \omega = 10$ MeV', linewidth=2)

	ax.legend(loc = 'best', frameon=False, fontsize = 18)

	fig.set_size_inches(8,8)

	ax.xaxis.set_tick_params(labelsize=20)
	ax.yaxis.set_tick_params(labelsize=20)

	plt.savefig('plots/' + plot_name + '.pdf', format='pdf')

def nmax_plot (nmax, energy):

	dir = os.path.abspath('..' + '/plots/')

	fig = plt.figure(frameon=False)

	ax = fig.add_subplot(111)

	ax.set_xlabel('$\hbar \omega$', fontsize = 22)
	ax.set_ylabel('Energy [MeV]', fontsize = 22)

	fig.suptitle('2 particles - trap $\hbar \omega$ = 10 MeV', fontsize = 22)

	plot_name = 'nmax' + str(nmax) + '_results'

	methods  = ['k-', 'r-',  'g-', 'b-', 'y-', 'm-', 'c-', 'r--', 'g--', 'b--']

	hw = [10, 15, 20, 25, 30, 35, 40, 45, 50, 60]

	energy2 = [27.76, 28.20, 30.08, 32.40, 35.35, 38.75, 42.45, 46.35, 50.40, 58.80]

	energy4 = [26.91, 26.75, 27.20, 28.42, 30.30, 32.33, 34.55, 36.96, 39.55, 45.09]

	plt.subplot(111).plot(hw, energy, methods[3], label = '${N}_{max} =$ ' + str(nmax), linewidth=2)

	ax.legend(loc = 'best', frameon=False, fontsize = 18)

	fig.set_size_inches(8,8)

	ax.xaxis.set_tick_params(labelsize=20)
	ax.yaxis.set_tick_params(labelsize=20)

	plt.savefig('plots/' + plot_name + '.pdf', format='pdf')


if __name__ == '__main__':

	energy2 = [27.76, 28.20, 30.08, 32.40, 35.35, 38.75, 42.45, 46.35, 50.40, 58.80]
	energy4 = [26.91, 26.75, 27.20, 28.42, 30.30, 32.33, 34.55, 36.96, 39.55, 45.09]

	onebody_plot2 ()
	onebody_plot4 ()
	nmax_plot (2, energy2)
	nmax_plot (4, energy4)
