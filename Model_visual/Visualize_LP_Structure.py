from scipy.sparse import coo_matrix
from pyomo.environ import *
from pyomo.util.model_size import build_model_size_report as model_size
import multiprocessing as mp
import time
import itertools
import numpy as np
def time_string(seconds):
	"""Returns time in seconds as a string formatted HHHH:MM:SS."""
	s = int(round(seconds))	 # round to nearest second
	h, s = divmod(s, 3600)	# get hours and remainder
	m, s = divmod(s, 60)  # split remainder into minutes and seconds
	return '%4i:%02i:%02i' % (h, m, s)


def extract_coeff(c,i, all_var):
	constr = c.split(" ")
	output = []
	for var in constr:
		if len(var) > 1:
			try:
				output.append((i,all_var[var]))

			except KeyError:
				pass
	return output
def get_struct(model,filename = None):
	all_var = {}
	xlab = []
	xloc = []
	index_val = 0
	start_overall = time.perf_counter()
	for c in model.component_data_iterindex(Var):
		if c[0][0] not in xlab:
			xlab.append(c[0][0])
			xloc.append(index_val)
		all_var[c[1].to_string()] = index_val
		# all_var.append(c[1].to_string())
		index_val += 1
	print("%s variables consolidated in %s" % (index_val, time_string(time.perf_counter() - start_overall)))

	start = time.perf_counter()
	row_data = []
	column_data = []
	input_val = []
	ylabels = []
	yloc = []
	num_constr = model_size(model).overall.constraints
	for c, i in zip(model.component_data_iterindex(Constraint), range(num_constr)):
		name = c[0][0]
		grouping = name.split("_")[0]
		if grouping not in ylabels:
			ylabels.append(grouping)
			yloc.append(i)
		constr = str(c[1]._body)
		constr = constr.split(" ")
		for var in constr:
			if len(var) > 1:
				try:
					row_data.append(all_var[var])
					column_data.append(i)
					input_val.append(1)
				except KeyError:
					pass
	print("%s constraints evaluated in %s" % (num_constr, time_string(time.perf_counter() - start)))
	Problem_Structure = coo_matrix((input_val,(column_data,row_data)),shape = (num_constr,len(all_var)))
	end_struct = time.perf_counter()
	print("Model structure generated in %s\n\n" % time_string(end_struct- start_overall))

		#def plot_structure(struct_matrix):
	import matplotlib.pyplot as plt
	a = plt.figure()
	plt.spy(Problem_Structure, markersize=.5, aspect='auto', color='red')
	if filename is None: a.savefig("LP_Structure.png")
	else: a.savefig(filename+".png")

	return Problem_Structure

