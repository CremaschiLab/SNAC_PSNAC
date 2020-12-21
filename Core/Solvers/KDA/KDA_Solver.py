import sys
import os
from pyomo.environ import *
from pyomo.opt import SolverFactory
import itertools
from pyutilib.misc import Options
import time as timer
import pdb
import gc
import re

class KDA_Solution:
	
	def __init__(self, knapsack_object, options, fixed_parameters={}):
		
		### Initialize Knapsack
		self. KDA_object =  knapsack_object
		
		### Initialize fixed parameters
		self.fixed_parameters = fixed_parameters	
		
		### Solve knapsack
		solution = self.kda_solver(options)
		
		self.output = solution
		
		
	def kda_solver(self, options):
		
		### Initialize Metrics
		problem_count = 0
		
		## Start Solution Timer
		start_time = timer.clock()

		## Set the number of time steps total
		Time_Periods = knapsack_object.Time_Steps
		
		## initialize time
		time_idx = 0
		
		### Iterate over time
		while time_idx < len(Time_Periods):
		
			### Generate list of this time period nodes
			ct_node = [node for node in self.knapsack_object.G.node if self.knapsack_object.G.node[node]['Time']=Time_Periods[time_idx]
			
			### For all sub-problems at t=t
			for nodes in ct_nodes:
				
				### Generate the Knapsack/EV problem
			
		