import sys
import os
import itertools
import copy
import math
import multiprocessing as mp
import Core.Solvers.MTSSP.M2S_item as M2S_item
from pyomo.environ import *
from pyomo.opt import SolverFactory
import Core.scenario_class as SC
from Core.Solvers.MSSP import defunction as MSSP
import pdb



def EOSS_PRDP_Solve(MD,FD):
	### Generate Scenarios
	NP = len(MD._data['product'][None])	
	NT = len(MD._data['trial'][None])
	
	Scenarios = itertools.product(range(NT + 1), repeat = NP)
	
	### Parallelization Parameters
	np = mp.cpu_count()
	
	gsize = math.ceil((NT+1)**NP/np)
		
	### Divide Objects into pool tasks
	setlist = _grouper(Scenarios,gsize)
	
	pool = mp.Pool(np)
	
	ENPV_Results = [pool.apply_async(Scen_Solve, args=(i,MD,FD)) for i in setlist]
	pool.close()
	pool.join()
	
	#ENPV_Results = 0
	#for i in setlist:
		#ENPV_Results = ENPV_Results + Scen_Solve(i,MD,FD)
	#ENPV = ENPV_Results	
	### Combine Results
	ENPV = 0
	
	for i in ENPV_Results:
		ENPV += i._value
		
	return ENPV

def Scen_Solve(S_List, MD,FD):
	ENPV_Set = 0
	### Build Initial Model
	model = build_model(MD)
	#print(FD)
	for s in S_List:
		
		if s != None:
			ss = tuple(s)
			#print(ss)
			### Definition
			product = MD._data['product'][None]
			trials = MD._data['trial'][None]
			
			## Scenario Specific Data
			sinf = SC.scenario(ss,MD._data['probability'],product,trials)
			
			Success = {}
			for i in product:
				model.Success[i] = sinf.success[product.index(i)]
			
			
			### Fix Fixed Variables
			for (i,j,t,dec) in FD:
				for rlzn in FD[(i,j,t,dec)]:
					if rlzn_check(ss,rlzn,product):
						#print(i,j,t+1,dec)
						model.Decision_X[i,j,t+1].value = dec
						model.Decision_X[i,j,t+1].fixed = True
			
			### Solve Model
			opt = SolverFactory("cplex")
			results = opt.solve(model)
			model.solutions.load_from(results)
			
			
			
			ENPV_Set = ENPV_Set + sinf.probability * float(results['Problem'][0]['Lower bound'])
			
			model.unfix_all_vars()
		
	return ENPV_Set

def rlzn_check(Outcome, rlzn, product):
	if len(rlzn) == 0:
		return True
	
	breal = 0
	for (i,j,real) in rlzn:
		if Outcome[product.index(i)] >= j and real == 1:
			breal += 1
		elif Outcome[product.index(i)] < j and real == 0:
			breal += 1
		else:
			return False
	
	if breal == len(rlzn):
		return True
	else:
		return False
	
def _grouper(ibl, n, fillvalue=None):
	args = [iter(ibl)]*n
	return itertools.zip_longest(fillvalue = fillvalue, *args)
	
def build_model(MD):
	from Core.Solvers.EOSS.SingleScenario_Model import SingleScenario as SS
	import Core.Solvers.MTSSP.M2S_item as M2S_item
	
	##Set product
	product = MD._data['product'][None]
		
	##Set stage_gate
	stage_gate = MD._data['trial'][None]
	
	## Set time step
	time_step = MD._data['time_step'][None]
	
	##Set resource type
	resource_type = MD._data['resource_type'][None]
	
	## Set duration
	duration = MD._data['trial_duration']
	
	## Set trial cost
	trial_cost = MD._data['trial_cost']
	
	## Set Discount Values
	gammaL={}
	for items in MD._data['gammaL']:
		gammaL[items[0]] = MD._data['gammaL'][items]
	
	gammaD = {}
	for items in MD._data['gammaD']:
		gammaD[items[0]] = MD._data['gammaD'][items]
	
	## Set Maximum Revenue
	revenue_max ={}
	for items in MD._data['maximum_revenue']:
		revenue_max[items[0]] = MD._data['maximum_revenue'][items]
		
	## Set Last Trial
	last_trial = len(stage_gate)
	
	last_time_step = len(time_step)
	
	## Calculate running rev
	rev_run = M2S_item.calc_rr(revenue_max,gammaL,duration, product, stage_gate, time_step)
		
	##Calculate open rev  
	rev_open = M2S_item.calc_openrev(revenue_max,gammaL,duration, product, stage_gate, time_step, last_time_step)
	
	##Calculate Discounting Factor
	discounting_factor = M2S_item.calc_discounting_factor(revenue_max,gammaL,trial_cost, product, stage_gate, last_time_step)
	
	resource_max = {}
	for items in MD._data['max_resource']:
		resource_max[items[0]] = MD._data['max_resource'][items]
		
	resource_required = {}
	resource_required = MD._data['resource_requirement']
	
	success = {}
	for i in product:
		success[i] = 0
	
	model = SS(product,stage_gate,time_step,resource_type,\
		resource_max,gammaL,gammaD,duration,\
		trial_cost,resource_required, revenue_max,\
		success,last_time_step, last_trial,\
		rev_run, rev_open, discounting_factor)
	
	return model
