import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from coopr.opt import SolverFactory, SolverManagerFactory
from pyutilib.misc import Options
import time
import gc

class warmstart_class:
	def __init__(self, model_data, kda_results, output_directory, mipgap=.001):
		### Process model_data for deterministic solve
		from MTSSP import PRDP_Data_Processing
		self.problem_data = PRDP_Data_Processing.MTSSP_PRDP_Data_Processing(model_data._data)
	
		self.results = kda_results.output['results']
		self.sp_realizations = kda_results.output['sub_problem_realizations']
		
		opt = SolverFactory("cplex")

		options = Options()
		opt.options.mip_tolerances_mipgap = mipgap
		
		### Generate Non-Anticipativity Constraints
		self.PRDP_NAC_Generation()
		self.Scen_Spec_Parms()
	
		### Generate Model
		from MSSP.defunction import de
		wrm_model_start_time = time.clock()
		model = de(self.problem_data.product,self.problem_data.stage_gate,self.problem_data.time_step,self.problem_data.resource_type,self.problem_data.SS,self.problem_data.resource_max, self.problem_data.gammaL,self.problem_data.gammaD,self.problem_data.duration,self.problem_data.trial_cost,self.problem_data.resource_required, self.problem_data.revenue_max,self.pb, self.problem_data.success,self.problem_data.Last_Time_Step, self.problem_data.last_trial, self.problem_data.running_revenue, self.problem_data.open_revenue, self.problem_data.discounting_factor, self.phi, self.phii, self.phij, self.outcome)
	
		### Create Instance
		wrm_instnce_strt_time = time.clock()
		instance = model.create()
		
		del model
		gc.collect()
		
		for s in self.problem_data.SS:
			### Calculate X
			xbox = self.Calculate_X(self.problem_data.List_of_Scenarios[s])
	
			### Fix X Values in Instance
			for i in self.problem_data.product:
				for j in self.problem_data.stage_gate:
					for t in self.problem_data.stage_gate:
						idx = self.problem_data.product.index(i)
						jdx = self.problem_data.stage_gate.index(j)
						tdx = self.problem_data.stage_gate.index(t)
						if xbox[idx][jdx][tdx]:
							instance.Decision_X[i,j,t,s].value == 1
		
		
	
	
		### Solve Model
		wrm_strt_time = time.clock()
		output = opt.solve(instance, warmstart=True)
		wrm_fnsh_time = time.clock()
		
		WS_Solve_Time = wrm_fnsh_time - wrm_strt_time
		WS_InstanceGen_Time =  wrm_strt_time - wrm_instnce_strt_time
		WS_ModelCreate_Time = wrm_instnce_strt_time - wrm_model_start_time
		Warmstart_Total_Time = wrm_fnsh_time - wrm_model_start_time
		
		### Clear Solution Variables
		del instance
		del output
		
		model_start_time = time.clock()
		
		### Recreate Model and Solve for Deterministic Time
		model = de(self.problem_data.product,self.problem_data.stage_gate,self.problem_data.time_step,self.problem_data.resource_type,self.problem_data.SS,self.problem_data.resource_max, self.problem_data.gammaL,self.problem_data.gammaD,self.problem_data.duration,self.problem_data.trial_cost,self.problem_data.resource_required, self.problem_data.revenue_max,self.pb, self.problem_data.success,self.problem_data.Last_Time_Step, self.problem_data.last_trial, self.problem_data.running_revenue, self.problem_data.open_revenue, self.problem_data.discounting_factor, self.phi, self.phii, self.phij, self.outcome)
		
		instnce_strt_time = time.clock()
		
		### Create Instance
		instance = model.create()
		
		### time and solve DE_Version
		strt_time = time.clock()
		de_results = opt.solve(instance)
		fnsh_time = time.clock()
		
		NWS_Solve_Time = fnsh_time - strt_time
		NWS_InstanceGen_Time = strt_time - instnce_strt_time
		NWS_ModelCreate_Time = instnce_strt_time - model_start_time
		Total_Time = fnsh_time - model_start_time
		
		### Load Results
		instance.load(de_results)
		
		### Transform Results
		transformed_results = instance.update_results(de_results)
		
		### Write File
		self.Output_Write(transformed_results, Warmstart_Total_Time,WS_Solve_Time,WS_InstanceGen_Time,WS_ModelCreate_Time,Total_Time, NWS_Solve_Time,NWS_InstanceGen_Time,NWS_ModelCreate_Time, output_directory)
		
		
	def PRDP_NAC_Generation(self):	
		#######################################################################
		### Generate Non-Anticipativity Constraints
		#######################################################################
		OC = {}
		for s in self.problem_data.SS:
			OC[s] = [] 
			for i in self.problem_data.product:
				OC[s].append(self.problem_data.List_of_Scenarios[s].outcome[self.problem_data.product.index(i)])
		
		self.phi= {}
		self.phii= {}
		self.phij ={}

		print "Generating Non-Anticipativity Constraints"
	
		for s in self.problem_data.SS:
			for sp in self.problem_data.SS:
				if sp > s:
					for i in self.problem_data.product:
						OCtest = list(OC[s])
						OCtest[self.problem_data.product.index(i)] += 1
						OCtest2 = list(OC[s])
						OCtest2[self.problem_data.product.index(i)] += -1
						if OCtest == OC[sp]:
							trl = OC[s][self.problem_data.product.index(i)] + 1
							self.phi[(s,sp)] = 1
							self.phii[(s,sp)] = i
							self.phij[(s,sp)] = trl
						if OCtest2 == OC[sp]:
							trl = OC[sp][self.problem_data.product.index(i)] + 1
							self.phi[(s,sp)] = 1
							self.phii[(s,sp)] = i
							self.phij[(s,sp)] = trl
		
	def Scen_Spec_Parms(self):
		self.pb = {}
		self.outcome = {}
			
		for s in self.problem_data.SS:
			self.pb[s] = self.problem_data.List_of_Scenarios[s].probability
			self.outcome[s] = self.problem_data.List_of_Scenarios[s].outcome
	
	def Calculate_X(self, scenario):
		 ### Generate box to store results
		xbox = []
		for i in self.problem_data.product:
			jbox = []
			for j in self.problem_data.stage_gate:
				tbox = [0] * len(self.problem_data.time_step)
				jbox.append(tbox)
			xbox.append(jbox)
		
		t = 0 
		
		while t < len(self.problem_data.time_step):	 					
			if t == 0:
		
				current_sp = '0'
				### The time zero decision is the same for all (root node sub-problem)
				for items in self.results[0]['0']:
					ii = self.problem_data.product.index(items[0])
					jj = self.problem_data.stage_gate.index(items[1])
					xbox[ii][jj][0] = 1
				t += 1
		
			else:
				### Try to match the current SP
				try: 
					self.results[t]
					if current_sp in self.results[t]:
						for items in self.results[t][current_sp]:
							ii = self.problem_data.product.index(items[0])
							jj = self.problem_data.stage_gate.index(items[1])
							xbox[ii][jj][t] = 1
						t += 1
						
					else:
					### If match use add decision
						searching = 1
						while searching == 1:
							for new_sp in self.results[t]:
								
								### IF no match, pop the entries (see if there was a realization)
								old_sp = new_sp.rsplit(".", 1)
								old_sp = old_sp[0]
								if current_sp == old_sp:
									matches = 0
							
									### If match, (should be min of two) see which sp_realization matches scenario outcome
									for selections in self.sp_realizations[new_sp]:
																			
										if scenario.outcome[selections[0]] > selections[1] and self.sp_realizations[new_sp][selections] == 1:
											matches += 1
										
										elif scenario.outcome[selections[0]] == selections[1] and self.sp_realizations[new_sp][selections] == 0:
											matches += 1
									
									
								
									if matches == len(self.sp_realizations[new_sp]):
										current_sp = new_sp
										for items in self.results[t][current_sp]:
											ii = self.problem_data.product.index(items[0])
											jj = self.problem_data.stage_gate.index(items[1])
											xbox[ii][jj][t] = 1
										t += 1
									
										searching = 0
									
														
								
				except:
					t += 1			
		return xbox

	def Output_Write(self,transformed_results, Warmstart_Total_Time,WS_Solve_Time,WS_InstanceGen_Time,WS_ModelCreate_Time,Total_Time, NWS_Solve_Time,NWS_InstanceGen_Time,NWS_ModelCreate_Time, output_directory):
		### Generate New File Name
		save_file =  "Warmstart" + "_" + "Output" 
		
		### Open save file
		if not os.path.exists(output_directory):
			os.makedirs(output_directory)
		
		f = open(os.path.join(output_directory, save_file),	"w")
		
		wrmstrt_tt = 'Warmstart Total Time:' + ' ' + str(Warmstart_Total_Time)
		f.write(wrmstrt_tt + '\n')
		
		wrmstrt_st = 'Warmstart Solve Time:' + ' ' + str(WS_Solve_Time)
		f.write(wrmstrt_st + '\n')
		
		wrmstrt_mt = 'Warmstart Model Time:' + ' ' + str(WS_ModelCreate_Time)
		f.write(wrmstrt_mt + '\n')
		
		wrmstrt_it = 'Warmstart Instance Time:' + ' ' + str(WS_InstanceGen_Time)
		f.write(wrmstrt_it + '\n')
		
		non_warmstart_time = 'Non-Warmstart Total Time: ' + str(Total_Time)
		f.write(non_warmstart_time + '\n')
		
		nwrmstrt_st = 'Non-Warmstart Solve Time:' + ' ' + str(NWS_Solve_Time)
		f.write(nwrmstrt_st + '\n')
		
		nwrmstrt_mt = 'Non-Warmstart Model Time:' + ' ' + str(NWS_ModelCreate_Time)
		f.write(nwrmstrt_mt + '\n')
		
		nwrmstrt_it = 'Non-Warmstart Instance Time:' + ' ' + str(NWS_InstanceGen_Time)
		f.write(nwrmstrt_it + '\n')
	
		de_results_write = str(transformed_results)
		f.write(de_results_write + '\n')
		
		f.close()
		
