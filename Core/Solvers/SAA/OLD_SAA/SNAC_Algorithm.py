import os
import sys

### Other Python Libraries
import itertools
import pdb
from operator import itemgetter
import multiprocessing as mp

### MST Imports
from .NAC_Graph import NAC_Graph
from .Kruskal_Class import Kruskal_MST as MST

###  imports
import Core.Solvers.SAA.SAA_NAC as SAA_NAC

def NAC_Generator(Uncertain_Parameters, Scenarios, Scenario_Outcomes):
	"""
		Uncertain Parameters - a tuple (i.e. ('A','B','C')) of the Uncertain Parameter Classes
								for the problem	(Uncertain Parameter class is in UP_Class.py)
		
		Scenarios-  a list of scenario numbers i.e. [1,2,3,...S]
		
		Scenario_Realizations is a dictionary mapping S to its outcome 
			i.e. { 1:(1,0,1)}, notice the realization of uncertain 
			parameters is a vector. The vector length should be equal to the 
			number of uncertain parameters
			
	""" 
	
	Parameters = [str(Uncertain_Parameters[i].Name) for i in range(len(Uncertain_Parameters))]	
	####################################################################
	### 				 Generate set of order rules 
	####################################################################
	
	Order_Rules = [[r,r.ROrder] for r in Uncertain_Parameters if r.ROrder != []]
	
	####################################################################
	### 				Generate set of possible cuts 
	####################################################################
	cuts = []
	
	"""
		This code loops over each uncertain parameter, if the uncertain
		parameter is realized instantaneously then the code appends the
		name of the parameter to the list of cuts. If the realization is
		gradual the code generates a set of cuts based on the way the 
		parameter is realized. Order Rules are updated to reflect the 
		gradual cuts.
	"""
	
	for up in Uncertain_Parameters:
		if up.Realization_Type == 'instant':
			cuts.append(up.Name)
		else:
			for grs in up.GROrderSets:
				cname = up.Name + '_' + str(grs)
				cbefore = [up.Name + '_' + str(pp) for pp in up.GROrderSets[:grs]]
				cuts.append(cname)
				Order_Rules.append([cname, cbefore])
	
	####################################################################
	### 					Generate initial MST  
	####################################################################
	
	### Initialize Graph Object
	graph = NAC_Graph()

	### Add Verticies to Graph
	for s in Scenarios:
		graph.add_vertex(s)
	
	### Generate All Edges
	graph.all_edges() 
	
	### Calculate MST using Kruskal Algorithm
	MST_Graph = MST(graph,Scenario_Outcomes)
	
	mst = MST_Graph.mst
	
	####################################################################
	### 					Generate cut sets  
	####################################################################
	C_List = itertools.permutations(cuts,len(cuts)-1)
	C_List = tuple(C_List)
	
	####################################################################
	### 					Generate NACs  
	####################################################################
	
	for c in C_List:
	
		### Check the validity of the cut list
		Valid = C_List_validation(c,Order_Rules)

		if Valid == True:
		
			sets = [Scenarios,]
			i = 0
	
			while i < len(c):
				
				### This is for instantaneous realizations
				if c[i] in Parameters:
					idx = Parameters.index(c[i])
					new_sets = []
					
					### Generate Multi-Process Object and Run MST Subset Gen
					if len(sets) < mp.cpu_count():
						np = len(sets)
					else:
						np = mp.cpu_count()
						
					### Create Pool
					pool = mp.Pool(np)
		
					### Run Pool
					results = [pool.apply_async(Instant_Parallel_MST, args=(ss,idx,Scenario_Outcomes,mst)) for ss in sets]
					pool.close()
					pool.join()
					
					### Consolodate Results
					mst = mst.union(*[results[i][0]._value for i in range(len(results))])	
									
					sets = [*results[i][1]._value for i in range(len(results))]
				
				### This is for gradual realizations
				else:
					cvar = c[i].split('_')
					if cvar[0] in Parameters:
						
						idx = Parameters.index(cvar[0])
						new_sets = []
						
						### Generate Multi-Process Object and Run MST Subset Gen
						if len(sets) < mp.cpu_count():
							np = len(sets)
						else:
							np = mp.cpu_count()
							
						### Create Pool
						pool = mp.Pool(np)
			
						### Run Pool
						results = [pool.apply_async(Gradual_Parallel_MST, args=(ss,idx,Scenario_Outcomes,mst)) for ss in sets]
						pool.close()
						pool.join()
						
						### Consolodate Results
						mst = mst.union(*[results[i][0]._value for i in range(len(results))])
							
						sets = [*results[i][1]._value for i in range(len(results))]
				i += 1
			
			
				
		
	pdb.set_trace()

def Instant_Parallel_MST(ss,idx,Scenario_Outcomes,mst):
	sortable = [(s,Scenario_Outcomes[s][idx]) for s in ss]
	subsets = groupby_func(sortable, key=itemgetter(1))
	new_sets.append(subsets)
			
	### Generate Multi-Process Object and Run MST Subset Gen
	
	if len(sets) < mp.cpu_count():
		np = len(sets)
	else:
		np = mp.cpu_count()
	### Create Pool
	pool = mp.Pool(np)

	### Run Pool
	results = [pool.apply_async(Subset_MST, args=(ss,Scenario_Outcomes,mst)) for ss in subsets]
	pool.close()
	pool.join()
	
	NAC_to_add = {}
	NAC_to_add = NAC_to_add.union(*[results[i]._value for i in range(len(results))])

	return NAC_to_add, subsets

def Gradual_Parallel_MST(ss,idx,Scenario_Outcomes,mst):
	subsets = []
	subset1 = [s for s in ss if str(scenario_outcomes[s][idx]) in list(cvar[1])]
	subset2 = [s for s in ss if str(scenario_outcomes[s][idx]) not in list(cvar[1])]
	if subset2 != []:
		new_sets.append(subset2)
		subsets.append(subset2)

	if len(list(cvar[1])) > 1:
		sortable = [(s,scenario_outcomes[s][idx]) for s in subset1]
		subset3 = groupby_func(sortable, key=itemgetter(1))
		if subset3 != []:
			new_sets.append(subset3)
			subsets.append(subset3)
	else:
		if subset1 != []:
			new_sets.append(subset1)
			subsets.append(subset1)
	
	### generate multi-process object and run mst subset gen
	
	if len(sets) < mp.cpu_count():
		np = len(sets)
	else:
		np = mp.cpu_count()
	
	### create pool
	pool = mp.pool(np)
	
	
	### run pool
	results = [pool.apply_async(subset_mst, args=(ss,scenario_outcomes,mst)) for ss in subsets]
	pool.close()
	pool.join()
	
	NAC_to_add = {}
	NAC_to_add = NAC_to_add.union(*[results[i]._value for i in range(len(results))])

	return NAC_to_add, subsets
	
def Subset_MST(subset,Scenario_Realizations,mst):
	existing_connection = []
	
	#### Get Old Edges
	for (s,sp) in mst:
		if (s,sp) in set(itertools.permutations(subset,2)):
			existing_connection.append((s,sp))
	
	if len(existing_connection) < len(subset) - 1:
		### Create Graph for subset1
		graph = NAC_Graph({})
		
		for s in subset:
			graph.add_vertex(s)
		graph.all_edges()
		
		### Find MST
		MST_subset = MST(graph, Scenario_Realizations, existing_connection)
		
		added_NACs = MST_subset.mst
	else:
		added_NACs = {}
		
	return added_NACs
	
def C_List_validation(C_List,Order_Rules):
	Valid = True
	truth_table = []
	for rules in Order_Rules:
		for i in rules[1]:
			if str(rules[0]) in C_List and str(i) in C_List and C_List.index(str(rules[0])) > C_List.index(str(i)):
				truth_table.append(True)
			elif str(rules[0]) not in C_List:
				truth_table.append(True)
			else:
				truth_table.append(False)
	if all(truth_table) and truth_table[0] == True:
		pass
	else:
		Valid = False
		return Valid
			
	return Valid

def groupby_func(data, key=itemgetter(0)):
	### Sort data
	data = sorted(data, key=key)
	
	### Generate subsets
	subset = [list(group) for k,group in itertools.groupby(data,key)]
	
	### Convert to useable form
	subset_return = []
	for subs in subset:
		subset_return.append([a[0] for a in subs])
		
	return subset_return
	
	 
