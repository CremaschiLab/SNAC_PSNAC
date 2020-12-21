import sys
import os
from coopr.opt import SolverFactory, SolverManagerFactory, SolverStatus, TerminationCondition, SolutionStatus
from pyutilib.misc import Options
from defunction import SingleScenario
import pdb
import gc


def EOSS_PRDP_Solve(s, problem_data, fixed_parameters,output_directory):
	opt = SolverFactory("cplex")

	options = Options()
	opt.options.mip_tolerances_mipgap = .000001
	opt.options.mip_tolerances_absmipgap = .000001
	
	SSsuccess = {}
	
	### Redefine Success
	for i in problem_data.product:
		SSsuccess[i] = problem_data.success[(i,s)]	
		
	### Redefine Outcome
	SSoutcome = problem_data.List_of_Scenarios[s].outcome
		
	### Redefine Probability
	SSProbability = problem_data.List_of_Scenarios[s].probability
		
	model = SingleScenario(problem_data.product,problem_data.stage_gate,problem_data.time_step,problem_data.resource_type,problem_data.resource_max,problem_data.gammaL,problem_data.gammaD,problem_data.duration,problem_data.trial_cost,problem_data.resource_required, problem_data.revenue_max,SSProbability, SSsuccess,problem_data.Last_Time_Step, problem_data.last_trial, problem_data.running_revenue, problem_data.open_revenue, problem_data.discounting_factor, SSoutcome)
		
	instance = model.create()
	
	###################################################################
	### Determine which fixed parameters are applicable to this problem
	###################################################################
	new_fixed_parameters = EOSS_Fixed_Parameters(fixed_parameters,SSoutcome)
	
	for items in new_fixed_parameters:	
		var_i = items[0]
		var_j = items[1]
		var_t = items[2] + 1
		decision = items[3]
		instance.Decision_X[var_i,var_j,var_t] = decision
		instance.Decision_X[var_i,var_j,var_t].fixed = True
		
		
	del model
	gc.collect()
	
	instance.preprocess()	
	results= opt.solve(instance)
	instance.load(results)	
	"""
	save_file = str(s) 
	
	### Open save file
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)
		
	f = open(os.path.join(output_directory, save_file),	"w")
	
	transformed_results = instance.update_results(results)
	tr = str(transformed_results)
	f.write(tr + '\n')
	f.close()
	"""	
	
	if results.solver.status == SolverStatus.ok and results.solver.termination_condition == TerminationCondition.optimal:
		return [s,results.solution.objective['__default_objective__']['Value']]
	else:
		save_file = "Scenario ", s
		if not os.path.exists(output_directory):
			os.makedirs(output_directory)
		
		f = open(os.path.join(output_directory, save_file),	"w")
	
		transformed_results = instance.update_results(results)
		tr = str(transformed_results)
		f.write(tr + '\n')
		f.close()
		exit()
		
	


def EOSS_Fixed_Parameters(fixed_parameters,SSoutcome):
	new_fixed_parameters = []
	for itms in fixed_parameters:
		for jtms in fixed_parameters[itms]:
			if len(jtms) == 0:
				new_fixed_parameters.append(itms)
			else:
				cntr = 0
				for ktms in jtms:
					if ktms[2] == 0:
						if SSoutcome[ktms[0]] == ktms[1]:
							cntr += 1
					else:
						if SSoutcome[ktms[0]] > ktms[1]:
							cntr += 1
				if cntr == len(jtms):
					new_fixed_parameters.append(itms)
	
	
	return new_fixed_parameters
