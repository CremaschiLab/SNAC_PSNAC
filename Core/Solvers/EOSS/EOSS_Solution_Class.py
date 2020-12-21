import sys
import os
from coopr.opt import SolverFactory, SolverManagerFactory
import itertools
from pyutilib.misc import Options
import time as timer
import pdb
import gc
import re
from Core.Solvers.MTSSP import MTSSP_Item
import multiprocessing as mp

class EOSS:
	def __init__(self,model_data,solver,solve_ouput,mipgap=.001,_opts=[], fixed_parameters = {}):
	
		self.mipgap = mipgap
		self.solver = solver
		self._opts = _opts
		
		### Generate Parameters
		problem_data = MTSSP_Item.MTSSP_Data_Processing(model_data)
		
		### Implement the solver
		
		
	def EOSS_Solver(self, problem_data, solve_output, _opts, fixed_parameters):
		
		### Set Solve set. This is always going to be the entire scenario set!!!
		if problem_data.model_type == 'PRDP':
			from EOSS_PRDP import EOSS_PRDP_Solve as function_name
			Scenario_Set = problem_data.SS
				
		################################################################
		###                  Solve Problems
		################################################################
		
		
		### Determine the number of cores
		process_count = mp.cpu_count()
		
		### Generate Pool
		pool = mp.Pool(processes = process_count)
		
		results = [pool.apply(function_name, args=(x,problem_data, fixed_parameters[x])) for x in Scenario_Set]
		
		
	
		
