import os
import sys
from Solver import solve_function
import Core.DataImport.parse_data_cmds
import Core.DataImport.import_data_class as import_data_class
import time
import itertools
import random
import math
import Core.scenario_class as scenario_class
from Core.Solvers.SAA.SNAC_Alg import NAC_Generator as NAC
from Core.Solvers.SAA.SNAC_Alg_original import NAC_Generator as NAC_original

from Core.Solvers.SAA.UP_Class import Uncertain_Parameter as UP
from Core.Solvers.MSSP.defunction import SAA as SAA
import Core.Solvers.MTSSP.M2S_item as M2S_item
from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
import statistics
import Core.Solvers.MSSP.defunction as defunction

from pyomo.core import Var
from pyomo.core import Constraint
from pyutilib.services import TempfileManager
import copy
import gc

import pdb
import json

def main():
	### Location and date information
	current_directory = os.path.dirname(os.path.realpath(__file__))
	current_date = time.strftime('%m_%d_%Y', time.gmtime())

	### Set Results File Output Directory
	output_directory = current_directory + '/Solutions/' + 'SAA_DE' + '_' + current_date + '/'

	### Loop Info
	# files = ['modeldata.dat','modeldatab.dat','modeldata3.dat']
	files = ['modeldata_4.dat','modeldata_10.dat','modeldata3.dat','modeldata3_5.dat','modeldata4_3.dat','modeldata4.dat','modeldata4_5.dat','modeldata5.dat','modeldata6.dat']
	# files = ['modeldata4_3.dat','modeldata4.dat','modeldata4_5.dat','modeldata5.dat','modeldata6.dat']
	# files = ['modeldata6.dat']
	# files = ['modeldata_4.dat']
	# Specify sample size
	# Number_of_Scenarios = {'modeldatab.dat': [.3, .5, .7, 1.0], 'modeldata.dat': [.3, .5, .7, 1.0],
						   # 'modeldata3.dat': [0.05, 0.10, 0.15, 0.20, 0.25, .3, .4, .5, 0.6, 0.8, 1.0]}
	# Number_of_Scenarios = {'modeldata6.dat':[.125, 0.25]}
	# Number_of_Scenarios = {'modeldata4.dat':[0.5]}
						   # 'modeldata_4.dat':[.75]
	# Number_of_Scenarios = {'modeldata_4.dat':[.75],'modeldata_10.dat':[.24],'modeldata3.dat':[.09375],'modeldata3_5.dat':[.192],'modeldata4_3.dat':[.148148148],'modeldata4.dat':[.5],'modeldata4_5.dat':[.0384],'modeldata5.dat':[.0625, 0.25, 0.5, 0.75, 1],'modeldata6.dat':[.125, 0.5, 0.75]}
	Number_of_Scenarios = {'modeldata_4.dat':[.75],'modeldata_10.dat':[.24],'modeldata3.dat':[.09375],'modeldata3_5.dat':[.192],'modeldata4_3.dat':[.148148148],'modeldata4.dat':[.5],'modeldata4_5.dat':[.0384],'modeldata5.dat':[.0625],'modeldata6.dat':[.125]}

	# Solve Count -- The specified number of time we solve the model.
	# scount = {'modeldata4.dat': 1}
	scount = {'modeldata.dat':1, 'modeldata_4.dat':1,'modeldata_10.dat':1,'modeldata3.dat':1,'modeldata3_5.dat':1,'modeldata4_3.dat':1,'modeldata4.dat':1,'modeldata4_5.dat':1,'modeldata5.dat':1,'modeldata6.dat':1}

	import_record_sceanrios = 'record_subset_scenarios_10_29_2020.json'
	with open(import_record_sceanrios) as f:
		Subset_Scenarios_record = dict([tuple((tuple(x[0]), x[1])) for x in json.loads(f.read())])

	for fl in files:
		print("Starting File " + str(fl))

		##### Import the data into the system
		model_data = get_datafile(fl)

		### Define number of products and number of trials
		NP = len(model_data._data['product'][None])

		NT = len(model_data._data['trial'][None])

		Product = model_data._data['product'][None]
		Trials = model_data._data['trial'][None]
		Probability = model_data._data['probability']

		### Generate all possible outcomes
		Outcomes = itertools.product(range(NT + 1), repeat=NP)
		Outcomes = list(Outcomes)

		### Generate Full Scenario Set
		scenario = 1
		Scenario_Objects = {}
		Scenario_Outcomes = {}
		Scenario_List = []

		for s in Outcomes:
			Scenario_Objects[scenario] = scenario_class.scenario(s, Probability, Product, Trials)
			Scenario_Outcomes[scenario] = s
			Scenario_List.append(scenario)
			scenario += 1

		Full_Sample_Size = len(Outcomes)
		Full_Sample_Size = math.ceil(Full_Sample_Size)
		

		n = 0
		nmax = 30
		sub_model_Max = 0
		sub_model_AVG = 0
		
		sub_model_Max_alter = 0
		sub_model_AVG_alter = 0
		
		full_model_Max = 0
		full_model_AVG = 0

		reducedfull_model_Max = 0
		reducedfull_model_AVG = 0
		NAC_AVG = 0
		NAC_MAX = 0
		
		NAC_AVG_alter = 0
		NAC_MAX_alter = 0
		
		while n < nmax:
			sub_model_NAC_count = 0
			sub_model_NAC_count_alter = 0
			full_model_NAC_count = 0
			reducedfull_model_NAC_count = 0

			for ns in Number_of_Scenarios[fl]:

				### Determine Sample Size
				Sample_Size = ns * len(Outcomes)
				Sample_Size = math.ceil(Sample_Size)

				### Select Scenarios
				Subset_Scenarios = Subset_Scenarios_record[fl,n,ns]
				print(len(Subset_Scenarios))

				### Normalize Probability
				total_probability = 0
				for s in Subset_Scenarios:
					total_probability += Scenario_Objects[s].probability

				S_Probability = {}
				for s in Subset_Scenarios:
					S_Probability[s] = Scenario_Objects[s].probability / total_probability

				Pb = {}
				for s in Scenario_List:
					Pb[s] = Scenario_Objects[s].probability

				### Generate NACs For Subset
				# print('Generating NACs for Scenarios' + str(Subset_Scenarios))
				Uncertain_Parameters = []
				for i in Product:
					Param = UP(str(i), 'endogenous', 'gradual', range(NT + 1), [], range(NT + 1))
					Uncertain_Parameters.append(Param)

				Uncertain_Parameters = tuple(Uncertain_Parameters)

				########################################################################################################
				################################### NAC GENERATION 1 ###################################################

				
				########################################################################################################
				################################### NAC GENERATION alternative ###################################################
				Ptime = time.process_time()
				Snac_StartTime = time.time()
				snacIter = 0
				while snacIter < scount[fl]:
					NAC_Set_alter = NAC(Uncertain_Parameters, Subset_Scenarios, Scenario_Outcomes)

					NAC_Set_alter = list(NAC_Set_alter)

					phI_IJ_alter = {}
					for (s, sp) in NAC_Set_alter:
						i = 0
						phI_IJ_alter[(s, sp)] = []
						while i < len(Product):
							phList = [Scenario_Outcomes[s][i], Scenario_Outcomes[sp][i]]
							if phList[0] != phList[1]:
								phList.sort()
								Differentiation = (Product[i], phList[0] + 1)
								phI_IJ_alter[(s, sp)].append(Differentiation)
							i += 1

					snacIter += 1
				Snac_Total_alter = time.time() - Snac_StartTime
				Snac_process_time = time.process_time() - Ptime

				## record NAC pairs
				NAC_AVG_alter += len(NAC_Set_alter)/1
				if len(NAC_Set_alter) > NAC_MAX_alter:
					NAC_MAX_alter = len(NAC_Set_alter)
					
				print('NAC_AVG_alter',NAC_AVG_alter,'NAC_MAX_alter',NAC_MAX_alter)
				list_phijj_alter = [x for x in phI_IJ_alter]
				print('list_phijj_alter',len(list_phijj_alter))

				########################################################################################################
				################################### NAC GENERATION 2 ###################################################


				########################################################################################################
				################################### NAC GENERATION 3 ###################################################
				### Generate NACs for Full Set

				####################################################################
				###                      Generate Model
				####################################################################
				Success = {}

				for s in Scenario_List:
					for i in Product:
						Success[(i, s)] = Scenario_Objects[s].success[Product.index(i)]

				GammaD = {}
				GammaL = {}
				revenue_max = {}

				resource_max = {}
				for items in model_data._data['max_resource']:
					resource_max[items[0]] = model_data._data['max_resource'][items]

				## Set Discount Values
				for items in model_data._data['gammaL']:
					GammaL[items[0]] = model_data._data['gammaL'][items]

				for items in model_data._data['gammaD']:
					GammaD[items[0]] = model_data._data['gammaD'][items]

				## Set Maximum Revenue
				for items in model_data._data['maximum_revenue']:
					revenue_max[items[0]] = model_data._data['maximum_revenue'][items]

				## Calculate running rev
				rev_run = M2S_item.calc_rr(revenue_max, GammaL, model_data._data['trial_duration'], Product, Trials,
										   model_data._data['time_step'][None])

				##Calculate open rev
				rev_open = M2S_item.calc_openrev(revenue_max, GammaL, model_data._data['trial_duration'], Product,
												 Trials, model_data._data['time_step'][None],
												 model_data._data['time_step'][None][-1])

				##Calculate Discounting Factor
				discounting_factor = M2S_item.calc_discounting_factor(revenue_max, GammaL,
																	  model_data._data['trial_cost'], Product, Trials,
																	  model_data._data['time_step'][None][-1])

				########################################################################################################
				################################## Create Model 1 ######################################################
				opt = SolverFactory("cplex")
				### generate number of NAC constriants
				
				########################################################################################################
				################################## Create Model 1 alternative ######################################################
				Ptime = time.process_time()
				Model1_StartTime = time.time()
				m1 = 0
				while m1 < scount[fl]:
					sub_model_alter = SAA(Product, Trials, model_data._data['time_step'][None],
									model_data._data['resource_type'][None], Subset_Scenarios, \
									resource_max, GammaL, GammaD, model_data._data['trial_duration'],
									model_data._data['trial_cost'], \
									model_data._data['resource_requirement'], revenue_max, S_Probability, \
									Success, len(model_data._data['time_step'][None]), Trials[-1], rev_run, rev_open,
									discounting_factor, \
									phI_IJ_alter, Scenario_Outcomes)
					m1 += 1
				Model1TotalTime_alter = time.time() - Model1_StartTime
				model_process_time = time.process_time() - Ptime
				
				##### Solve Sub Model Model
				Ptime = time.process_time()
				Model1Solve_StartTime = time.time()
				s1 = 0
				while s1 < scount[fl]:
					sub_results = opt.solve(sub_model_alter)
					s1 += 1
				Model1Solve_TotalTime_alter = time.time() - Model1Solve_StartTime
				model_solve_process_time = time.process_time() - Ptime
				sub_model_alter.solutions.load_from(sub_results)
				
				
				### generate number of NAC constriants
				sub_model_NAC_count_alter = len(sub_model_alter.NAC_Constraint) + len(sub_model_alter.NAC2_Constraint) + len(sub_model_alter.NAC3_Constraint)
				print('sub_model_NAC_count_alter',sub_model_NAC_count_alter)
				
				sub_model_AVG_alter += sub_model_NAC_count_alter/1
				
				if sub_model_NAC_count_alter > sub_model_Max_alter:
					sub_model_Max_alter = sub_model_NAC_count_alter

				print('sub_model_Max_alter',sub_model_Max_alter,'sub_model_AVG_alter',sub_model_AVG_alter)
				### generate number of NAC constriants
				
				
				########################################################################################################
				################################## Create Model 2 ######################################################


				########################################################################################################
				########################################################################################################

				### Generate New File Name
				save_file = 'SNAC_Alter_' + str(fl) + "_" + str(ns) + "_" + str(nmax) + "_iterations"

				### Open save file
				if not os.path.exists(output_directory):
					os.makedirs(output_directory)
					
					
				if n == 0:
					f = open(os.path.join(output_directory, save_file), "w")
					Header = "%-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s" %('SAA_par_ENPV',
					'SNAC_par_Time', 'SNAC_par_Model_Time', 'SNAC_par_Time(P)', 'SNAC_par_Model_Time(P)', 
					'SNAC_par_sol_Time', 'SNAC_par_sol_Time(P)', 
					'sub_model_par_AVG','sub_model_par_Max',
					'NAC_Count_AVG_par', 'NAC_Count_MAX_par')
					f.write(Header)
					f.write('\n')
				else:
					f = open(os.path.join(output_directory, save_file), "a")
				

				
				### Generate file contents
				RESULTS = "%-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s %-30s" %(str(round(sub_model_alter.Expected_NPV(), 4)),
					str(round(Snac_Total_alter, 4)),  str(round(Model1TotalTime_alter, 4)), str(round(Snac_process_time, 4)), str(round(model_process_time, 4)), 
					str(round(Model1Solve_TotalTime_alter, 4)), str(model_solve_process_time),
					str(round(sub_model_AVG_alter, 4)),str(round(sub_model_Max_alter, 4)),
					str(NAC_AVG_alter), str(NAC_MAX_alter)
					)

				f.write(RESULTS)
				f.write('\n')

			n += 1
			
		f.write('\n')
		f.write('----------' + import_record_sceanrios + '---------------')
		f.write('\n')
		f.close()

def get_datafile(modelfile):
	### Find the file to import
	problem_file_directory = os.path.dirname(os.path.realpath(__file__)) + '/Problem Files/'

	### If the file is in the problem files folder then the filename is the sub-directory plus the filename
	if os.path.isfile(os.path.join(problem_file_directory, modelfile)):
		import_file_name = os.path.join(problem_file_directory, modelfile)

	### Otherwise assume the file is in the current directory
	else:
		import_file_name = modelfile

	### Import command for file
	file_reader = ['import', import_file_name]

	### Import Data
	print("Importing data from " + str(import_file_name))
	model_data = import_data_class.Data_Collection(file_reader)

	return model_data


if __name__ == "__main__":
	main()
