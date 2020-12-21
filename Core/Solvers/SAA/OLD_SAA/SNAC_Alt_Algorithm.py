import os
import sys

### Other Python Libraries
import itertools
import pdb
from operator import itemgetter
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
import copy
import time
import math
import pickle
import zlib

### MST Imports
from .NAC_Graph import NAC_Graph
from .Kruskal_Class import Kruskal_MST as MST
from .CSO_Class import CSO 

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
			g = 0
			while g < len(up.GROrderSets)-1:
				cname = up.Name + '_' + str(up.GROrderSets[g]) + '_' + str(up.GROrderSets[g + 1])
				cbefore = [up.Name + '_' + str(pp) + '_' + str(pp+1) for pp in up.GROrderSets[:g]]
				cuts.append(cname)
				Order_Rules.append([cname, cbefore])
				g += 1
	
	####################################################################
	### 					Generate cut sets  
	####################################################################
	utype = [up.Realization_Type for up in Uncertain_Parameters]
	if all(utype[i] == 'gradual' for i in range(len(utype))):
		Rlists = []
		for up in Uncertain_Parameters:
			k = 0
			while k < len(up.GROrderSets) - 1:
				Rlists.append(up.Name)
				k += 1
		C_List = {tuple(p) for p in next_permutation(Rlists)}
		
		
		### Convert C_List to useable

		### create pool
		pool = mp.Pool(np)
		
		gsize = math.ceil(math.factorial(len(cuts))/np)
		
		### Divide Objects into pool tasks
		setlist = _grouper(C_List,gsize)
		
		Valid_C_List = [pool.apply_async(Conv_C_List, args=(i,Uncertain_Parameters)) for i in setlist]
		pool.close()
		pool.join()
		
	else:
		C_List = itertools.permutations(cuts,len(cuts))
	
		### create pool
		pool = mp.Pool(np)
		
		gsize = math.ceil(math.factorial(len(cuts))/np)
		
		### Divide Objects into pool tasks
		setlist = _grouper(C_List,gsize)
		del C_List
		
		Valid_C_List = [pool.apply_async(Val_Par, args=(i,Order_Rules)) for i in setlist]
		pool.close()
		pool.join()
	
	####################################################################
	###     Create Files to Hold Large Data Bits
	####################################################################
	n_count = 0
	C_List_Compressed = {}
	C_List_Compressed[0] = []
	for i in Valid_C_List:
		f = pickle.loads(zlib.decompress(i._value))
		if len(f) != 0:
			name =  'CList_Hold_' + str(n_count)
			CComp = CSO(name,f, Scenarios)
			C_List_Compressed[0].append(zlib.compress(pickle.dumps(CComp)))
			n_count += 1
		
	####################################################################
	### 					Generate NACs  
	####################################################################
	k = 0
	mst = {}
	current_date = time.strftime('%m_%d_%Y', time.gmtime())
	output_directory = current_directory + '/Solutions/' + 'NAC_Generation' + '/' + str(len(Parameters))+'_Parameters' + '_' + current_date + '/'
	
	
	while k < len(cuts)-1:
		mst_hold = {}
		for m in C_List_Compressed[k]:
			obj = pickle.loads(zlib.decompress(m))
			### create pool
			pool = mp.Pool(np)
			
			### Divide Objects into pool tasks
			setlist = _grouper(obj.CLists,np)
			setlist = list(setlist)
			
			results = [pool.apply_async(Unnest_Par_Func, args=(obj.sets, Scenario_Outcomes, Parameters, mst,i, k)) for i in setlist]
			pool.close()
			pool.join()
			
			#for i in setlist:
				#rresults, rresults2 = Unnest_Par_Func(obj.sets,Scenario_Outcomes,Parameters,mst,i,k)
			
			### Update Results
			mst = set(mst)


			nmst = mst.union(*[set(results[i]._value[1]) for i in range(len(results))])

			mst_hold[C_List_Compressed[k].index(m)] = nmst	
			sets = dict()
			for i in range(len(results)):
				for j in results[i]._value[0]:
					sets[j] = results[i]._value[0][j]
			
			obj.sets = sets
			try:
				C_List_Compressed[k+1].append(zlib.compress(pickle.dumps(obj))) 
			except:
				C_List_Compressed[k+1] = []
				C_List_Compressed[k+1].append(zlib.compress(pickle.dumps(obj))) 
			
		new_mst = set()
		nwt = new_mst.union(*[mst_hold[fm] for fm in range(len(C_List_Compressed[k]))])
		mst = nwt
		k += 1
	
	### Check for mst of full set
	
	sets = [s for s in Scenarios]
	NACs = Subset_MST(sets,Scenario_Outcomes,mst)
	NACs = set(NACs).union(mst)
	return NACs	
				
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

def Unnest_Par_Func(sets,SR,Parameters,mst,Valid_C_List,k):
	new_sets = {}
	mst = set(mst)
	for v in Valid_C_List:
		if v != None:
			new_sets[tuple(v)] = Combine_Group(sets[tuple(v)],SR,list(v[k+1:]), Parameters)
			print(sets[tuple(v)])
			print(list(v[k+1:]))
			for ss in new_sets[tuple(v)]:
				NAC_add = Subset_MST(ss,SR,mst)
				mst = mst.union(set(NAC_add))
	return new_sets, mst

def Parallel_Function(sets, SR, Parameters, mst, Valid_C_List,k):
	new_sets = {}
	mst = set(mst)
	for v in Valid_C_List:
		if v != None:	
			new_sets[v] = Combine_Group(sets[v],SR,list(v[k+1:]), Parameters)
		
			setlist = tuple(_grouper(new_sets[v],nsets))
			

			np = 3
			nsets = int(math.ceil(len(new_sets[v])/np))
			
			setlist = tuple(_grouper(new_sets[v],nsets))
					
			### create pool
			pool = mp.Pool(len(setlist))
					
			### run pool	
			results = [pool.apply_async(Subset_MST, args=(ss,SR,mst)) for ss in setlist]
			pool.close()
			pool.join()
			NAC_to_add = set()
			NAC_to_add = NAC_to_add.union(*[set(results[i]._value) for i in range(len(results))])
			
			mst = mst.union(set(NAC_to_add))

	return  new_sets, mst

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

def cmp(a,b):
	 return (a > b) - (a < b)
	 
def reverse(seq, start, end):
        # seq = seq[:start] + reversed(seq[start:end]) + \
        #       seq[end:]
        end -= 1
        if end <= start:
            return
        while True:
            seq[start], seq[end] = seq[end], seq[start]
            if start == end or start+1 == end:
                return
            start += 1
            end -= 1
            
### I didn't write this code. I found it on stack overflow. 
### Hoping it works better than the set reduce function	
def next_permutation(seq, pred = cmp):
    """Like C++ std::next_permutation() but implemented as
    generator. Yields copies of seq."""
    
    if not seq:
        raise StopIteration
    try:
        seq[0]
    except TypeError:
        raise TypeError("seq must allow random access.")
    first = 0
    last = len(seq)
    seq = seq[:]
    # Yield input sequence as the STL version is often
    # used inside do {} while.
    yield seq[:]
    if last == 1:
        raise StopIteration
    while True:
        next = last - 1
        while True:
            # Step 1.
            next1 = next
            next -= 1
            if pred(seq[next], seq[next1]) < 0:
                # Step 2.
                mid = last - 1
                while not (pred(seq[next], seq[mid]) < 0):
                    mid -= 1
                seq[next], seq[mid] = seq[mid], seq[next]
                # Step 3.
                reverse(seq, next1, last)
                # Change to yield references to get rid of
                # (at worst) |seq|! copy operations.
                yield seq[:]
                break
            if next == first:
                raise StopIteration
    raise StopIteration 
	
def C_List_validation(C_List,Order_Rules):
	Valid = True
	truth_table = []
	
	if len(Order_Rules) == 0:
		return Valid
	
	### For all rules
	for rules in Order_Rules:
		for i in rules[1]:
			
			### Both cuts are in the cutlist, and the cut order holds
			if str(rules[0]) in C_List and str(i) in C_List and C_List.index(str(rules[0])) < C_List.index(str(i)):
				truth_table.append(True)
				
			### The cut that must come before is the first cut
			elif str(i) not in C_List:
				truth_table.append(True)
			else:
				truth_table.append(False)
	
	if all(truth_table) and truth_table[0] == True:
		pass
	else:
		Valid = False
		return Valid
			
	return Valid

def Val_Par(process_list,Order_Rules):
	rlist = []
	for i in process_list:
		if i != None:
			if C_List_validation(i,Order_Rules) == True:
				rlist.append(i)
	return zlib.compress(pickle.dumps(rlist))

def Conv_C_List(C_List, UP):
	UP_Names = [UP[i].Name for i in range(len(UP))]
	Return_Lists = []
	for i in list(C_List):
		if i != None:
			cnt = [1 for k in UP]
			temp_list = []
			for s in i:
				upidx = UP_Names.index(s)
				cname = str(UP[upidx].Name) + '_' + str(UP[upidx].GROrderSets[-(cnt[upidx] +1)])  + '_' + str(UP[upidx].GROrderSets[-cnt[upidx]])
				temp_list.append(cname)
				cnt[upidx] += 1
			Return_Lists.append(temp_list)	
	return zlib.compress(pickle.dumps(Return_Lists))				
	
def _grouper(ibl, n, fillvalue=None):
	args = [iter(ibl)]*n
	return itertools.zip_longest(fillvalue = fillvalue, *args)
