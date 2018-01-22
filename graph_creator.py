import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from itertools import product
from scipy.interpolate import interp2d
import numpy as np
import os

def onebody_plot ():

	dir = os.path.abspath('..' + '/plots/')

	fig = plt.figure(frameon=False)

	ax = fig.add_subplot(111)

	ax.set_xlabel('${N}_{max}$', fontsize = 22)
	ax.set_ylabel('Energy [MeV]', fontsize = 22)

	fig.suptitle('4 particles - trap $\hbar \omega$ = 10 MeV', fontsize = 22)

	plot_name = 'onebody_results'

	methods  = ['k-', 'r-',  'g-', 'b-', 'y-', 'm-', 'c-', 'r--', 'g--', 'b--']

	nmax = [2, 4, 6, 8]

	hw10 = [80.00, 80.00, 80.00, 80.00]
	hw15 = [84.38, 80.49, 80.00, 80.00]
	hw20 = [93.98, 83.32, 80.76, 80.15]
	hw25 = [106.2, 88.49, 82.84, 80.91]
	hw30 = [120.0, 95.34, 86.30, 82.57]
	hw35 = [134.7, 103.3, 90.86, 85.12]
	hw40 = [149.9, 112.1, 96.25, 88.46]
	hw45 = [165.6, 121.4, 102.3, 92.42]
	hw50 = [181.5, 131.2, 108.8, 96.89]
	hw60 = [214.0, 151.5, 122.8, 106.9]

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


if __name__ == '__main__':

	onebody_plot ()
