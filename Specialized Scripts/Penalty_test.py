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
	output_directory = current_directory + '/Solutions/' + 'Penalty_Test' + '_' + current_date + '/'
	files = ['modeldata.dat','modeldata3.dat','modeldata4.dat','modeldata5.dat','modeldata6.dat']
	
	for fl in files:
	
		##### Import the data into the system
		model_file = fl
		model_data = get_datafile(model_file)
		
		### Set output directory		
		ps_output_directory = output_directory +  '/' + str(fl) + '/'

		### Get base case KDA solutions
		kda_output_directory = ps_output_directory + '/' + 'KDA' + '/' + 'Original' + '/' 
		OKDA = KDA_Solve_Function(model_data, kda_output_directory, [], ['max_solve'])
		
		### Increments of 1 to 100:
		penval = range(1,501)
		
		for p in penval:
			
			### Solve KDA Penalty with penval
			nkda_output_directory = ps_output_directory + '/' + 'KDA' + '/' + str(p)
			NKDA = KDA_Solve_Function(model_data,nkda_output_directory, {}, ['max_solve','penalty'], penval = p)
			
			### Compare KDA Solutions
			comp_results = SolComp(OKDA,NKDA)
			
			### Write Results
			ENPV_Percent = (NKDA.ENPV - OKDA.ENPV)/ OKDA.ENPV
			filewriterdata = [fl,p,NKDA.ENPV,ENPV_Percent,comp_results,NKDA.Problem_Count, NKDA.Solve_Time]
			_write_file(filewriterdata,output_directory)
			
	pdb.set_trace()
			
			
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
	save_file =  "Penalty" + "_" + "Output" 
		
	if not os.path.isfile(os.path.join(directory,save_file)):
		
		### Open save file
		f = open(os.path.join(directory, save_file),	"w")
		
		headings = "%-25s %-25s %-25s %-25s %-25s %-25s %-25s" % ('Parent File', 'Penalty Value', 'ENPV', 'Percent Difference', 'Solution Change?','Problem Count', 'Solve Time')
		line_break = '-' * 25 * 10
		zero_line = "%-25s %-25s %-25s %-25s %-25s %-25s %-25s" % (data[0], str(data[1]), str(round(data[2],4)),str(round(data[3],4)),str(data[4]),str(data[5]),str(round(data[6],4)))
		
		f.write(headings + "\n")
		f.write(line_break + "\n")
		f.write(zero_line + "\n")
		f.close()
	else:
		f = open(os.path.join(directory, save_file),	"a")
		zero_line = "%-25s %-25s %-25s %-25s %-25s %-25s %-25s" % (data[0], str(data[1]), str(round(data[2],4)),str(round(data[3],4)),str(data[4]),str(data[5]),str(round(data[6],4)))
		f.write(zero_line + "\n")
		f.close()
		
def KDA_Solve_Function(model_data, output_directory, fixed_parameters={},_opts=[], solver='cplex', mipgap=.001, penval = 0):
	## Run KDA to get lower bound, this returns the result and the ENPV location
	import Core.Solvers.KDA.KDA_Solution_Class as Solve
	import Core.Solvers.KDA.Evaluate_KDA_PRDP as Evaluate_KDA
	
	Solution = Solve.KDA(model_data, solver, mipgap, output_directory,_opts,fixed_parameters, penval)
	
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
	
if __name__ == "__main__":
	main_func()
