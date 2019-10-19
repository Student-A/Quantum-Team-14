
from subprocess import DEVNULL, STDOUT, check_call
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.stats import rv_discrete
import time


def normalize(v):
    norm = np.linalg.norm(v, ord=1)
    if norm == 0:
        norm = np.finfo(v.dtype).eps
    return v / norm

def get_point(dim, steps):
	dim = 'N ' + str(dim) + '\n\n'
	step_details = 'steps ' + str(steps) + '\n\n'
	parse_file = open("parse_file.txt","w")

	with open("remainder.txt","r") as remainder:
		parse_file.write(dim)
		parse_file.write(step_details)
		parse_file.write(remainder.read())
		parse_file.close()
		remainder.close()
		check_call('./main', stdout=DEVNULL, stderr=STDOUT)

		# os.system('./main')
		time.sleep(1)
		probability_data = open('./Output/Probability_Distribution', 'r')
		xs = []
		ys = []
		zs = []

		for line in probability_data:
			data = line.split(' ')
			xs.append(float(data[0]))
			ys.append(float(data[1]))
			zs.append(float(data[2]))

		zipped = list(zip(xs, ys))
		point = np.random.choice(range(len(zipped)), 1, p=normalize(zs))[0]
		return zipped[point]


def show_distribution_plot():
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	probability_data = open('./Output/Probability_Distribution', 'r')

	xs = []
	ys = []
	zs = []

	for line in probability_data:
		data = line.split(' ')
		xs.append(float(data[0]))
		ys.append(float(data[1]))
		zs.append(float(data[2]))
	xs = np.array(xs)
	ys = np.array(ys)
	zs = np.array(zs)
	dz = zs.flatten()
	bottom = np.zeros_like(xs)
	ax.bar3d(xs, ys, bottom, 0.5,  0.5, dz, shade=True)
	plt.show()

print(get_point(5, 500))
show_distribution_plot()
