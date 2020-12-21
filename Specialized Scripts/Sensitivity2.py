import os
import sys
import time as timer
from Solver import solve_function
import Core.DataImport.parse_data_cmds 
import Core.DataImport.import_data_class as import_data_class
from Core.Bounding.Node_Decision_Tree import Decision_Tree as DT
import time
import itertools
import copy
import gc
import resource
import pdb

def get_datafile(modelfile):
	### Find the file to import
	problem_file_directory = os.path.dirname(os.path.realpath(__file__)) + '/Problem Files/'
		
	### If the file is in the problem files folder then the filename is the sub-directory plus the filename
	if os.path.isfile(os.path.join(problem_file_directory, modelfile)):
		import_file_name = os.path.join(problem_file_directory,modelfile)
			
	### Otherwise assume the file is in the current directory
	else:	
		import_file_name = modelfile
	
	### Import command for file	
	file_reader = ['import', import_file_name]
	
	### Import Data
	print("Importing data from " + str(import_file_name))    
	model_data = import_data_class.Data_Collection(file_reader)
	
	
	return model_data


def main_func():
	### Location and date information
	current_directory = os.path.dirname(os.path.realpath(__file__))
	current_date = time.strftime('%m_%d_%Y', time.gmtime())
	
	### Set Results File Output Directory
	output_directory = current_directory + '/Solutions/' + 'Sensitivity2' + '_' + current_date + '/'
	files = ['modeldata5.dat']
	parms = ['nresources', 'ntrials', 'nplanning_horizon','cost','revenue','penalty1','penalty2']
	
	parmmod = [.75,.9,1.1,1.25]
	PHmod = [ -1 ,1,2,3]
	Tmod = [1,2,3]
	Rmod = [1,2,3]
	
	### Extend formulation info (Number of Trials)
	TC_Trial_Mod = [[1.25,1.4,1.5],[1.15,1.2,1.6],[1.1,1.3,1.8],[1.2,1.35,1.5],[1.3,1.35,1.5],[1.15,1.4,1.6]]
	Res1_Trial_Mod = [[3,3,3],[2,3,3],[2,2,3],[3,3,3],[3,2,3],[2,3,3]]   
	Res2_Trial_Mod = [[2,2,3],[2,3,3],[3,2,3],[2,3,3],[3,3,3],[3,3,3]]
	Dur_Trial_Mod = [[3,4,4],[3,4,5],[4,4,5],[3,3,5],[4,4,5],[4,3,5]]
	Prob_Trial_Mod = [[.8,.85,.9],[.75,.8,.9],[.8,.8,.85],[.85,.9,.95],[.75,.9,.9],[.8,.85,.85]]
	
	### Extend formulation info (Number of Resources)
	Resource_Mimic_Resource_Mod = [[[2,1,1],[1,1,1],[1,2,1]],[[1,2,2],[2,1,2],[1,1,2]],[[2,2,1],[1,2,1],[1,1,1]],[[2,2,2],[1,1,2],[2,1,2]],[[2,1,1],[2,2,1],[1,2,1]],[[1,2,2],[2,2,2],[1,1,2]]]
	Resource_Change_Resource_Mod = [[[1,-1,-1],[0,0,1],[-1,1,1]],[[0,-1,1],[1,-1,1],[-1,0,0]],[[-1,-1,0],[1,0,-1],[1,1,0]],[[0,1,0],[-1,-1,1],[1,1,0]],[[0,1,1],[-1,-1,-1],[1,-1,-1]],[[0,0,-1],[1,1,1],[1,-1,0]]]
	Max_Resource_Mod = [3,4,3]
	
	for fl in files:
	
		##### Import the data into the system
		model_file = fl
		model_data = get_datafile(model_file)
		
		### Set output directory		
		ps_output_directory = output_directory +  '/' + str(fl) + '/'

		### Get base case KDA solutions
		kda_output_directory = ps_output_directory + '/' + 'KDA' + '/' + 'Original' + '/' 
		OKDA = KDA_Solve_Function(model_data, kda_output_directory, [], ['max_solve'])
				
		for p in parms:
			
			
			### Modify problem data based on parm
			if p == 'cost':
				
				### Loop over different variations
				for m in parmmod:
					
					### Make copy of modeldata
					mmodeldata = copy.deepcopy(model_data)
					
					### Update change parameter	
					for (i,j) in mmodeldata._data['trial_cost']:
							mmodeldata._data['trial_cost'][(i,j)] = m * mmodeldata._data['trial_cost'][(i,j)]
		
					### Solve problem
					newKDAOutput = ps_output_directory + '/' + 'KDA' + '/' + p + '/' + 'mod' + str(m)
					
					if not os.path.exists(newKDAOutput):
						os.makedirs(newKDAOutput)
						
					ModKDA = KDA_Solve_Function(mmodeldata, newKDAOutput, [], [])
					
					### Compare KDA Solutions
					comp_results = SolComp(OKDA,ModKDA)
					
					### Solve Deterministic Problem
					mssp_output_directory = ps_output_directory + '/' + 'MSSP' + '/' + p + '/' + 'mod' + str(m)
					MSSP_Results = deterministic_solve_function(mmodeldata, mssp_output_directory)
					
					
					### Write Results
					ENPV_Percent = (ModKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
					filewriterdata = [fl,p,m,ModKDA.ENPV,ENPV_Percent,comp_results,ModKDA.Problem_Count, ModKDA.Solve_Time, MSSP_Results.ENPV, MSSP_Results.Total_Time]
					_write_file(filewriterdata,output_directory)
			
			elif p == 'revenue':
				
				### Loop over different variations
				for m in parmmod:
					
					### Make copy of modeldata
					mmodeldata = copy.deepcopy(model_data)
					
					### Update change parameter	
					for i in mmodeldata._data['maximum_revenue']:
							mmodeldata._data['maximum_revenue'][i] = m * mmodeldata._data['maximum_revenue'][i]
		
					### Solve problem
					newKDAOutput = ps_output_directory + '/' + 'KDA' + '/' + p + '/' + 'mod' + str(m)
					
					if not os.path.exists(newKDAOutput):
						os.makedirs(newKDAOutput)
						
					ModKDA = KDA_Solve_Function(mmodeldata, newKDAOutput, [], [])
					
					### Compare KDA Solutions
					comp_results = SolComp(OKDA,ModKDA)
					
					### Solve Deterministic Problem
					mssp_output_directory = ps_output_directory + '/' + 'MSSP' + '/' + p + '/' + 'mod' + str(m)
					MSSP_Results = deterministic_solve_function(mmodeldata, mssp_output_directory)
					
					
					### Write Results
					ENPV_Percent = (ModKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
					filewriterdata = [fl,p,m,ModKDA.ENPV,ENPV_Percent,comp_results,ModKDA.Problem_Count, ModKDA.Solve_Time, MSSP_Results.ENPV, MSSP_Results.Total_Time]
					_write_file(filewriterdata,output_directory)
					
			elif p == 'penalty1':
				
				### Loop over different variations
				for m in parmmod:
					
					### Make copy of modeldata
					mmodeldata = copy.deepcopy(model_data)
					
					### Update change parameter	
					for i in mmodeldata._data['gammaL']:
							mmodeldata._data['gammaL'][i] = m * mmodeldata._data['gammaL'][i]
		
					### Solve problem
					newKDAOutput = ps_output_directory + '/' + 'KDA' + '/' + p + '/' + 'mod' + str(m)
					
					if not os.path.exists(newKDAOutput):
						os.makedirs(newKDAOutput)
						
					ModKDA = KDA_Solve_Function(mmodeldata, newKDAOutput, [], [])
					
					### Compare KDA Solutions
					comp_results = SolComp(OKDA,ModKDA)
					
					### Solve Deterministic Problem
					mssp_output_directory = ps_output_directory + '/' + 'MSSP' + '/' + p + '/' + 'mod' + str(m)
					MSSP_Results = deterministic_solve_function(mmodeldata, mssp_output_directory)
					
					
					### Write Results
					ENPV_Percent = (ModKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
					filewriterdata = [fl,p,m,ModKDA.ENPV,ENPV_Percent,comp_results,ModKDA.Problem_Count, ModKDA.Solve_Time, MSSP_Results.ENPV, MSSP_Results.Total_Time]
					_write_file(filewriterdata,output_directory)
					
			elif p == 'penalty2':
				
				### Loop over different variations
				for m in parmmod:
					
					### Make copy of modeldata
					mmodeldata = copy.deepcopy(model_data)
					
					### Update change parameter	
					for i in mmodeldata._data['gammaD']:
							mmodeldata._data['gammaD'][i] = m * mmodeldata._data['gammaD'][i]
		
					### Solve problem
					newKDAOutput = ps_output_directory + '/' + 'KDA' + '/' + p + '/' + 'mod' + str(m)
					
					if not os.path.exists(newKDAOutput):
						os.makedirs(newKDAOutput)
						
					ModKDA = KDA_Solve_Function(mmodeldata, newKDAOutput, [], [])
					
					### Compare KDA Solutions
					comp_results = SolComp(OKDA,ModKDA)
					
					### Solve Deterministic Problem
					mssp_output_directory = ps_output_directory + '/' + 'MSSP' + '/' + p + '/' + 'mod' + str(m)
					MSSP_Results = deterministic_solve_function(mmodeldata, mssp_output_directory)
					
					
					### Write Results
					ENPV_Percent = (ModKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
					filewriterdata = [fl,p,m,ModKDA.ENPV,ENPV_Percent,comp_results,ModKDA.Problem_Count, ModKDA.Solve_Time, MSSP_Results.ENPV, MSSP_Results.Total_Time]
					_write_file(filewriterdata,output_directory)
			
			elif p == 'nplanning_horizon':
				
				### Loop over the different variations
				for m in PHmod:
					
					### Make copy of modeldata
					mmodeldata = copy.deepcopy(model_data)
					
					### Update parameter
					OPH = len(mmodeldata._data['time_step'][None])
					NewPH = OPH + m
					mmodeldata._data['time_step'][None] = [i+1 for i in range(NewPH)]
					
					### Solve problem
					newKDAOutput = ps_output_directory + '/' + 'KDA' + '/' + p + '/' + 'mod' + str(m)
					
					if not os.path.exists(newKDAOutput):
						os.makedirs(newKDAOutput)
						
					ModKDA = KDA_Solve_Function(mmodeldata, newKDAOutput, [], [])
					
					### Compare KDA Solutions
					comp_results = SolComp(OKDA,ModKDA)
					
					### Solve Deterministic Problem
					mssp_output_directory = ps_output_directory + '/' + 'MSSP' + '/' + p + '/' + 'mod' + str(m)
					MSSP_Results = deterministic_solve_function(mmodeldata, mssp_output_directory)
					
					
					### Write Results
					ENPV_Percent = (ModKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
					filewriterdata = [fl,p,m,ModKDA.ENPV,ENPV_Percent,comp_results,ModKDA.Problem_Count, ModKDA.Solve_Time, MSSP_Results.ENPV, MSSP_Results.Total_Time]
					_write_file(filewriterdata,output_directory)
					
			elif p == 'ntrials':
				
				### Loop over variations
				for m in Tmod:		
					
					### Make copy of modeldata
					mmodeldata = copy.deepcopy(model_data)
					
					### Update parameter
					OT = len(mmodeldata._data['trial'][None])
					NewT = OT + m
					mmodeldata._data['trial'][None] = [i+1 for i in range(NewT)]
					
					### Update Other Parameters Modified By Parameter Modification
					for i in mmodeldata._data['product'][None]:
						for j in mmodeldata._data['trial'][None][OT:]:
							mmodeldata._data['trial_cost'][(i,j)] = TC_Trial_Mod[mmodeldata._data['product'][None].index(i)][mmodeldata._data['trial'][None][OT:].index(j)] * mmodeldata._data['trial_cost'][(i,j-1)]
							mmodeldata._data['trial_duration'][(i,j)] = Dur_Trial_Mod[mmodeldata._data['product'][None].index(i)][mmodeldata._data['trial'][None][OT:].index(j)] 
							mmodeldata._data['probability'][(i,j)] = Prob_Trial_Mod[mmodeldata._data['product'][None].index(i)][mmodeldata._data['trial'][None][OT:].index(j)]
							mmodeldata._data['resource_requirement'][(i,j,'Type1')] = Res1_Trial_Mod[mmodeldata._data['product'][None].index(i)][mmodeldata._data['trial'][None][OT:].index(j)]
							mmodeldata._data['resource_requirement'][(i,j,'Type2')] = Res2_Trial_Mod[mmodeldata._data['product'][None].index(i)][mmodeldata._data['trial'][None][OT:].index(j)]
					
					### Solve problem
					newKDAOutput = ps_output_directory + '/' + 'KDA' + '/' + p + '/' + 'mod' + str(m)
					
					if not os.path.exists(newKDAOutput):
						os.makedirs(newKDAOutput)
						
					ModKDA = KDA_Solve_Function(mmodeldata, newKDAOutput, [], [])
					
					### Compare KDA Solutions
					comp_results = SolComp(OKDA,ModKDA)
					
					### Solve Deterministic Problem
					mssp_output_directory = ps_output_directory + '/' + 'MSSP' + '/' + p + '/' + 'mod' + str(m)
					MSSP_Results = deterministic_solve_function(mmodeldata, mssp_output_directory)
					
					
					### Write Results
					ENPV_Percent = (ModKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
					filewriterdata = [fl,p,m,ModKDA.ENPV,ENPV_Percent,comp_results,ModKDA.Problem_Count, ModKDA.Solve_Time, MSSP_Results.ENPV, MSSP_Results.Total_Time]
					_write_file(filewriterdata,output_directory)
					
			else:
				
				### Loop over variations
				for m in Rmod:
					
					### Make copy of modeldata
					mmodeldata = copy.deepcopy(model_data)
					
					### Update Parameter
					OldResNum = len(mmodeldata._data['resource_type'])
					cnt = OldResNum + 1
					added_types = []
					
					while cnt <= OldResNum + m:
						added_types.append('Type' + str(cnt + 1))
						cnt += 1
					
					mmodeldata._data['resource_type'][None] += added_types
					
					### Update Parameters impacted by change
					for rtypes in added_types:
						for i in mmodeldata._data['product'][None]:
							for j in mmodeldata._data['trial'][None]:
								mimic = Resource_Mimic_Resource_Mod[mmodeldata._data['product'][None].index(i)][mmodeldata._data['trial'][None].index(j)][added_types.index(rtypes)]
								mimic = 'Type' + str(mimic)
								rtypemod = Resource_Change_Resource_Mod[mmodeldata._data['product'][None].index(i)][mmodeldata._data['trial'][None].index(j)][added_types.index(rtypes)]
								mmodeldata._data['resource_requirement'][(i,j,rtypes)] = mmodeldata._data['resource_requirement'][(i,j,mimic)] + rtypemod
						
						mmodeldata._data['max_resource'][(rtypes,)] = Max_Resource_Mod[added_types.index(rtypes)] 		
					
					
					### Solve problem
					newKDAOutput = ps_output_directory + '/' + 'KDA' + '/' + p + '/' + 'mod' + str(m)
					
					if not os.path.exists(newKDAOutput):
						os.makedirs(newKDAOutput)
						
					ModKDA = KDA_Solve_Function(mmodeldata, newKDAOutput, [], [])
					
					### Compare KDA Solutions
					comp_results = SolComp(OKDA,ModKDA)
					
					### Solve Deterministic Problem
					mssp_output_directory = ps_output_directory + '/' + 'MSSP' + '/' + p + '/' + 'mod' + str(m)
					MSSP_Results = deterministic_solve_function(mmodeldata, mssp_output_directory)
					
					
					### Write Results
					ENPV_Percent = (ModKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
					filewriterdata = [fl,p,m,ModKDA.ENPV,ENPV_Percent,comp_results,ModKDA.Problem_Count, ModKDA.Solve_Time, MSSP_Results.ENPV, MSSP_Results.Total_Time]
					_write_file(filewriterdata,output_directory)
		
					
						
						
					
					
def SolComp(OldSol, NewSol):
	### Default Value of Solution matching is True
	sol_match = True
	
	### For all the times in the original solution
	for t in OldSol.Results:
		
		### Is decision time is in the new solution
		if t in NewSol.Results:
			
			### Check to see if the sub-problems and decisions match
			for sp in OldSol.Results[t]:
				
				if sp != '0':
					
					### Find the realization for the selected sub-problem
					sprel = OldSol.Realizations[sp]
					
					### Determine the sub-problem in the new solution that matches
					spmatch = None
					for sps in NewSol.Realizations:
						if NewSol.Realizations[sps] == sprel:
							spmatch = sps
					
					if spmatch == None:
						sol_match = False
						return sol_match
					
					if spmatch in NewSol.Results[t]:
						if set(OldSol.Results[t][sp]) != set(NewSol.Results[t][spmatch]):
							sol_match = False
							return sol_match
					else:
						sol_match = False
						return sol_match
				
		else:
			sol_match = False
			return sol_match			
					
					
	return sol_match
	

def _write_file(data,directory):
	### Generate New File Name
	save_file =  "Sensitivity" + "_" + "Output" 
		
	if not os.path.isfile(os.path.join(directory,save_file)):
		
		### Open save file
		f = open(os.path.join(directory, save_file),	"w")
		
		headings = "%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s" % ('Parent File', 'Modification', 'Mod Value', 'ENPV', 'Percent Difference', 'Solution Change?','Problem Count', 'Solve Time', 'MSSP ENPV', 'MSSP Total Time')
		line_break = '-' * 25 * 10
		zero_line = "%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s" % (data[0], data[1], str(data[2]),str(round(data[3],4)),str(round(data[4],4)),str(data[5]),str(data[6]),str(round(data[7],4)), str(round(data[8],4)), str(round(data[9],4)))
		
		f.write(headings + "\n")
		f.write(line_break + "\n")
		f.write(zero_line + "\n")
		f.close()
	else:
		f = open(os.path.join(directory, save_file),	"a")
		zero_line = "%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s" % (data[0], data[1], str(data[2]),str(round(data[3],4)),str(round(data[4],4)),str(data[5]),str(data[6]),str(round(data[7],4)), str(round(data[8],4)), str(round(data[9],4)))
		f.write(zero_line + "\n")
		f.close()

def KDA_Solve_Function(model_data, output_directory, fixed_parameters={},_opts=[], solver='cplex', mipgap=.001 ):
	## Run KDA to get lower bound, this returns the result and the ENPV location
	import Core.Solvers.KDA.KDA_Solution_Class as Solve
	import Core.Solvers.KDA.Evaluate_KDA_PRDP as Evaluate_KDA
	
	Solution = Solve.KDA(model_data, solver, mipgap, output_directory,_opts,fixed_parameters)
	
	#### Calculate the equivalent ENPV
	results = Solution.output['results']
	sp_realizations = Solution.output['sub_problem_realizations']
	Evaluated_Solution = Evaluate_KDA.KDA_PRDP_results(model_data,results,sp_realizations, output_directory)
	
	### Write results to Consolidated file
	import Core.output_write as output_write
	output_write._write(Evaluated_Solution,Solution,output_directory,'kda', '')
	
	import Core.Bounding_Class as Bounding_Class
	return_value =  Bounding_Class.Bkda(Evaluated_Solution,Solution)
	
	return return_value
	
def deterministic_solve_function(model_data, output_directory, mipgap=.05):
	
	import Core.Solvers.MSSP.Deterministic_Solver as Solve
	MSSP_Solution = Solve.deterministic_PRDP_solve_with_return(mipgap, model_data, output_directory)
	
	return MSSP_Solution

if __name__ == "__main__":
	main_func()
