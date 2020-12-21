import os
import sys

### Other Python Libraries
import itertools
import pdb
from operator import itemgetter
import multiprocessing as mp
import copy
import time
import math
import pickle
import zlib
from collections import Counter as mset

### MST Imports
from .NAC_Graph import NAC_Graph
from .Kruskal_Class import Kruskal_MST as MST


def NAC_Generator(Uncertain_Parameters, Scenarios, Scenario_Outcomes, opts=[], current_directory=''):
	"""Uncertain Parameters - a tuple (i.e. ('A','B','C')) of the Uncertain Parameter Classes
								for the problem	(Uncertain Parameter class is in UP_Class.py)
		
		Scenarios-  a list of scenario numbers i.e. [1,2,3,...S]
		
		Scenario_Realizations is a dictionary mapping S to its outcome 
			i.e. { 1:(1,0,1)}, notice the realization of uncertain 
			parameters is a vector. The vector length should be equal to the 
			number of uncertain parameters
			
	""" 
	np = mp.cpu_count()
	Parameters = [str(Uncertain_Parameters[i].Name) for i in range(len(Uncertain_Parameters))]	
	####################################################################
	### 				 Generate set of order rules 
	####################################################################
	
	Order_Rules = [[r,r.ROrder] for r in Uncertain_Parameters if r.ROrder != []]

	###################################################################
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
			g = 0
			while g < len(up.GROrderSets)-1:
				cname = up.Name + '_' + str(up.GROrderSets[g]) + '_' + str(up.GROrderSets[g + 1])
				cbefore = [up.Name + '_' + str(pp) + '_' + str(pp+1) for pp in up.GROrderSets[:g]]
				cuts.append(up.Name)
				Order_Rules.append([cname, cbefore])
				g += 1
	
	
	####################################################################
	### 					Generate NACs  
	####################################################################
	k = 1
	mst = set()
	current_date = time.strftime('%m_%d_%Y', time.gmtime())
	output_directory = current_directory + '/Solutions/' + 'NAC_Generation' + '/' + str(len(Parameters))+'_Parameters' + '_' + current_date + '/'
	
	
	sets = {}
	sets[()] = [[j,] for j in Scenarios]
	
	if len(Scenarios) > 10000000:
		use_sets = True
	else:
		use_sets = False
		

	while k < len(cuts):
		
		if len(Scenarios) > 500:
			print('SNAC Iteration: ' + str(k))
		
		### Generate combination of length k from cut list 
		clists = set(itertools.combinations(cuts,k))
		clists = list(clists)
		# print('K',k,'clists',clists)
		### Get Added NACs for each item in the clists
		################################################################
		
		### create pool
		# pool = mp.Pool(16, maxtasksperchild=1000)
		counttt = mp.cpu_count()
		pool = mp.Pool(counttt)
		
		nsets = math.ceil(len(clists)/np)	
		
		### Divide Objects into pool tasks
		setlist = _grouper(list(clists),nsets)
		setlist = list(setlist)

		results = [pool.apply_async(clist_manager, args=(Scenario_Outcomes, Scenarios,Uncertain_Parameters, Parameters, mst,cuts,i,sets,use_sets)) for i in setlist]
		pool.close()
		pool.join()
		# results = {}
		#pdb.set_trace()
		
		# for i in setlist:
			# results[i], ww = clist_manager(Scenario_Outcomes, Scenarios, Uncertain_Parameters, Parameters, mst, cuts, i,sets,use_sets)
			
		# print(results)
		# print(results)
		# for i in clists:
			# results[i] = clist_manager(Scenario_Outcomes, Scenarios, Uncertain_Parameters, Parameters, mst, cuts, i,sets,use_sets)
		
		### Update Results
		# print(results)
		mst = set(mst)
		# for i in range(len(results)):
			# print(set(results[i]))

		try:

			# nmst = mst.union(*[set(results[i]._value[0]) for i in setlist])
			nmst = mst.union(*[set(results[i]._value[0]) for i in range(len(results))])
			# nmst = mst.union(*[set(results[i]) for i in setlist])
		except:
			pdb.set_trace()
			
		mst = nmst
		
		### Update Sets
		sets = {}
		for i in results:
			sets.update(i._value[1])


		k += 1
		
	sets = [s for s in Scenarios]
	NACs = Subset_MST(sets,Scenario_Outcomes,mst)
	NACs = set(NACs).union(mst)
	return NACs	
		
def clist_manager(Scenario_Outcomes, Scenarios, Uncertain_Parameters, Parameters, mst, cuts, i, sets, use_sets):
	new_set_hold = {}
	### Combine groups
	for v in i:
		if v != None: 
			if not use_sets:
				new_sets = [[j,] for j in Scenarios]
				Remaining_Realizations = copy.deepcopy(cuts)
				
				### This builds sets from the ground up, not ideal 
				for c in v:
					Remaining_Realizations.remove(c)
					
					## Convert to Group Combine Form (1_2_3, 1_1_2, 1_0_1)
					RR = copy.deepcopy(Remaining_Realizations)
					for up in Uncertain_Parameters:
						ll = 0
						if up.Realization_Type == 'gradual':
							for k in Remaining_Realizations:
								if k == str(up.Name):
									idx = RR.index(k)
									RR[idx] = RR[idx] + '_' + str(ll) + '_' + str(ll + 1)
									ll += 1	
				
					new_sets = Combine_Group(new_sets,Scenario_Outcomes,RR, Parameters)
			else:
				Remaining_Realizations = copy.deepcopy(cuts)
				
				for c in v:
					Remaining_Realizations.remove(c)
				
				RR = copy.deepcopy(Remaining_Realizations)
				for up in Uncertain_Parameters:
						ll = 0
						if up.Realization_Type == 'gradual':
							for k in Remaining_Realizations:
								if k == str(up.Name):
									idx = RR.index(k)
									RR[idx] = RR[idx] + '_' + str(ll) + '_' + str(ll + 1)
									ll += 1	
				
				for sl in sets:
					
					if len(list((mset(v) & mset(sl)).elements())) == len(v) - 1:
						## Convert this
						new_sets = Combine_Group(sets[sl],Scenario_Outcomes,RR, Parameters)
						break
					
				if 'new_sets' not in locals():
					print(mset(v))
				
			new_set_hold[tuple(v)] = copy.deepcopy(new_sets)
			for ss in new_sets:
				NAC_add = Subset_MST(ss,Scenario_Outcomes,mst)
				mst = mst.union(set(NAC_add))
				
		# print('ss',new_sets)
				
	return mst,new_set_hold
			
def Combine_Group(groups, SR, Remaining_Realizations, Parameters):
	return_groups = []
	dynamic_groups = copy.deepcopy(groups)
	
	
	for g in dynamic_groups:
		new_group = copy.deepcopy(g)
		for gg in dynamic_groups:
			if groups.index(g) < groups.index(gg):
				
				Match = Group_Compare(SR[g[0]],SR[gg[0]],Remaining_Realizations, Parameters)
				
				if Match: 
					new_group += copy.deepcopy(gg)
					dynamic_groups.remove(gg)
					
		
		### if we've added groups together
		return_groups.append(new_group)
				
	return return_groups

def Group_Compare(A,B,IP,UP):
	Match = True
	for condition in IP:
		if condition in UP:
			""" This implies that the uncertainty is realized instantaneously"""
			idx = UP.index(condition)
			if A[idx] != B[idx]:
				""" The condition is not met"""
				Match = False
				return Match
		else:
			""" These are gradual realizations of uncertainty"""
			cmap = condition.split('_')
			idx = UP.index(cmap[0])
			if str(A[idx]) <= cmap[1] and str(B[idx]) >= cmap[2]:
				Match = False
				return Match
			elif str(A[idx]) >= cmap[2] and str(B[idx]) <= cmap[1]:
				Match = False
				return Match
	return Match
	
def Subset_MST(subset,Scenario_Realizations,mst):
	list_added_NAC = set()
	
	existing_connection = []
		
	#### Get Old Edges
	for (v,vp) in mst:
		if (v,vp) in set(itertools.permutations(subset,2)):
			existing_connection.append((v,vp))
				
		
	#if len(existing_connection) < len(subset) - 1:
	### Create Graph for subset1
	graph = NAC_Graph({})
		
	for ss in subset:
		graph.add_vertex(ss)
	graph.all_edges()
		
	### Find MST
	MST_subset = MST(graph, Scenario_Realizations, existing_connection)
		
	added_NACs = MST_subset.mst
	#else:
		#added_NACs = {}
	
	list_added_NAC = list_added_NAC.union(set(added_NACs))
	
		
	return list_added_NAC
		
def _grouper(ibl, n, fillvalue=None):
	args = [iter(ibl)]*n
	return itertools.zip_longest(fillvalue = fillvalue, *args)
