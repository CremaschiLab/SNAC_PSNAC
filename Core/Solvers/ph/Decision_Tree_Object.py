import os
import sys
import networkx
import pdb
import matplotlib
import matplotlib.pyplot as plt
import itertools
import copy

class Multistage_PYSP_Object:
	def __init__(self,Fixed_Variables, ModelData):
		
		if Fixed_Variables == {('Drug2', 1, 0, 1): [set()], ('Drug3', 1, 0, 0): [set()], ('Drug3', 1, 1, 1): [set()], ('Drug1', 1, 0, 0): [set()], ('Drug2', 2, 2, 1): [{('Drug2', 1, 1)}]}:
			pdb.set_trace()
			
		## Convert Decisions into decision tree
		###############################################################
		## Convert Fixed Variables to Useable Format
		
		self.debug = False
				
		tree_map, tree_realizations = self.fixed_variable_conversion(Fixed_Variables)
		
		## Define Graph and Attributes
		
		self.G = networkx.DiGraph()
		
		self.build_scenario_network(tree_map,tree_realizations,ModelData)	
			
		self.Construct_PS(ModelData)
		
		self.add_scenarios(ModelData)
		
		
		#if Fixed_Variables == {('Drug1', 1, 0, 1):[set()], ('Drug2', 1, 0, 0):[set()], ('Drug2', 1, 1, 1):[set()], ('Drug1', 2, 2, 1):[{('Drug1', 1, 1)}],('Drug2', 2, 3, 0):[{('Drug2', 1, 1), ('Drug1', 1, 0)}]}:
			#from Core.Solvers.PH.Network_Plot import hierarchy_pos as hierarchy_pos
			#fig1 = plt.figure()
			#pos = hierarchy_pos(self.G,'R')
			#networkx.draw_networkx(self.G, pos=pos,with_labels = True)
			#print(tree_map)
			#print(tree_realizations)
			#plt.show()
			
		#from Core.Solvers.PH.Network_Plot import hierarchy_pos as hierarchy_pos
		#fig1 = plt.figure()
		#pos = hierarchy_pos(self.G,'R')
		#networkx.draw_networkx(self.G, pos=pos,with_labels = True)
		#print(tree_map)
		#print(tree_realizations)
		#plt.show()
		
	def build_scenario_network(self,tree_map,tree_realizations,MD):
		time_periods =[]
		self.stage_names = []
		nt = 0
		tmax = max([t for t in tree_map])

		for t in MD._data['time_step'][None]:
			
			if t <= tmax:
				if "T" + str(t) not in self.stage_names:
					self.stage_names.append("T" + str(t))
			
				if t == MD._data['time_step'][None][0]:
					
					## Add root node and add information
					if len(tree_map[t]) > 1:
						print("There is something Fishy Going On!")
						pdb.set_trace()
						
					for key in tree_map[t]:
						self.G.add_node("R")
						self.G.node['R']['Realizations'] = []
						self.G.node['R']['Fixed Variables'] = tree_map[t][key]
						self.G.node['R']['Active Periods'] = MD._data['time_step'][None][0]
						
											
				else:
					
					### There are no new nodes just nodes that persist
					previous_time_nodes = [n for n in self.G if self.G.node[n]['Active Periods'] == t-1]
					
					for nodes in previous_time_nodes:
						### Check there are realizations but not enough resources to start new trials
						
						## Determine the path between root node and previous time node
						rpath = networkx.all_simple_paths(self.G,'R',nodes)
						rpath = list(rpath)
							
						### This identifies the case where there is more than one path
						if len(rpath) > 1: 
							print('Error in Tree')
							exit()
						
						if nodes == 'R':
							cummulative_decisions = self.G.node['R']['Fixed Variables']
						else:
							all_decs = [self.G.node[nint]['Fixed Variables'] for nint in rpath[0]]
							cummulative_decisions = []
							for itr in all_decs:
								cummulative_decisions += itr
								
						realizations = []		
						for (ip,jp,tp,decs) in cummulative_decisions:
							if MD._data['trial_duration'][(ip,jp)] + tp == t and decs == 1:
								realizations.append([ip,jp])
						
						### If there aren't any realizations then the nodes persist
						if len(realizations) == 0:
							self.G.add_node('u' + str(nt))
							self.G.node['u' + str(nt)]['Realizations'] = self.G.node[nodes]['Realizations']
							### Check tree_map to see if there are fixed variables
							if t in tree_map:
								match = False
								for branch in tree_map[t]:
									if set(tree_realizations[branch]) == set(self.G.node[nodes]['Realizations']):	
										match = True
										self.G.node['u' + str(nt)]['Fixed Variables'] = tree_map[t][branch]
								if match == False:
									self.G.node['u' + str(nt)]['Fixed Variables'] = []
							else:		
								self.G.node['u' + str(nt)]['Fixed Variables'] = []
							self.G.node['u' + str(nt)]['Active Periods'] = t
							self.G.add_edge(nodes,'u' + str(nt), probability = 1)
						
							nt += 1
							
						### Otherwise we need to construct the approapriate number of nodes
						else:
							nnn = 2 ** len(realizations)
							otcs = itertools.product(range(2), repeat= len(realizations))
							otcs = list(otcs)
							add_list = []
							for itr in otcs:
								add_item = [(realizations[i][0],realizations[i][1],itr[i]) for i in range(len(itr))]
								add_list.append(add_item)
							
							node_add = 0
							while node_add < nnn:
				
								### add ps node
								self.G.add_node('u' + str(nt))
								self.G.node['u' + str(nt)]['Active Periods'] = copy.deepcopy(self.G.node[nodes]['Active Periods']) + 1
								
								### Determine realizations
								reals = set(add_list[node_add]).union(set(self.G.node[nodes]['Realizations']))
								self.G.node['u' + str(nt)]['Realizations'] = reals
								
								### Fixed Variables
								if t in tree_map:
									match = False
									for branch in tree_map[t]:
										if set(tree_realizations[branch]) == set(self.G.node['u' + str(nt)]['Realizations']):	
											match = True
											self.G.node['u' + str(nt)]['Fixed Variables'] = tree_map[t][branch]
									if match == False:
										self.G.node['u' + str(nt)]['Fixed Variables'] = []
								else:		
									self.G.node['u' + str(nt)]['Fixed Variables'] = []
									
								### attach to parent
								node_pb = self.calc_p(MD,set(add_list[node_add]))	
								self.G.add_edge(nodes,'u' + str(nt), probability=node_pb) 
								nt += 1
								node_add += 1
									
	def calc_p(self, MD, rlzn):
		pb = 1
		for (i,j,r) in rlzn:
			if r == 0:
				pb = pb * (1- MD._data['probability'][(i,j)])
			else:
				pb = pb * MD._data['probability'][(i,j)]
		
		return pb		
			
	def fixed_variable_conversion(self,Fixed_Variables):
		tree = {}
		tree_realizations = {}
		tmap = 0
		for (i,j,t,d) in Fixed_Variables:
			if t + 1 not in tree:
				tree[t+1] = {}
				for nd in range(len(Fixed_Variables[(i,j,t,d)])):
					tree[t+1][tmap] = []
					tree[t+1][tmap].append((i,j,t+1,d))
					
					tree_realizations[tmap] = Fixed_Variables[(i,j,t,d)][nd]
					tmap += 1
				
			else:
				for nd in range(len(Fixed_Variables[(i,j,t,d)])):
					appended = False
					for branch in tree[t+1]:
						if set(Fixed_Variables[(i,j,t,d)][nd]) == set(tree_realizations[branch]):
							tree[t+1][branch].append((i,j,t+1,d))
							appended = True
					if appended == False:
						tree[t+1][tmap] = []
						tree[t+1][tmap].append((i,j,t+1,d))
							
						tree_realizations[tmap] = Fixed_Variables[(i,j,t,d)][nd]
						tmap += 1
				
				
		
		return tree, tree_realizations
		
	def Construct_PS(self,ModelData):
		### Construct nodes that represent realizations of uncertainty that have no decision on them
		nc = 0
		node_list = [node for node in self.G.node if self.G.out_degree(node) != 0]
		
		for node in node_list:
				
			## Check to make sure the probability out is equal to one
			calc_p = 0
			for node_out in self.G.edge[node]:
				calc_p += self.G.edge[node][node_out]['probability']
				
			#If the probability isn't equal to one then determine how many ps nodes to add	
			if calc_p != 1:
				
				### Determine the difference in realization between node and node_out
				differences = []
				for node_out in self.G.edge[node]:
					diff = set(self.G.node[node_out]['Realizations']).difference(set(self.G.node[node]['Realizations']))
					diff = list(diff)
					differences.append(diff)
				
				### Calculate the number of new nodes based on the differences
				number_of_new_nodes = 2 ** len(differences[0]) - len(differences)
				otcs = itertools.product(range(2), repeat= len(differences[0]))
				otcs = list(otcs)
				
				ijlist = []
				for (i,j,r) in differences[0]:
					ijlist.append([i,j])
					
				### Generate all realizations and remove the ones that already have nodes
				add_list = []
				
				mod_diff = [set(differences[i]) for i in range(len(differences))]
				for itr in otcs:
					add_item = [(ijlist[i][0],ijlist[i][1],itr[i]) for i in range(len(itr))]
					if set(add_item) not in mod_diff:
						add_list.append(add_item)
				
			
				node_add = 0
				while node_add < number_of_new_nodes:
					
					### add ps node
					self.G.add_node('up' + str(nc))
					self.G.node['up' + str(nc)]['Fixed Variables'] = []
					self.G.node['up' + str(nc)]['Active Periods'] = copy.deepcopy(self.G.node[node]['Active Periods']) + 1
					
					### Determine realizations
					reals = set(add_list[node_add]).union(set(self.G.node[node]['Realizations']))
					self.G.node['up' + str(nc)]['Realizations'] = reals
					
					### attach to parent
					node_pb = self.calc_p(ModelData,set(add_list[node_add]))	
					self.G.add_edge(node,'up' + str(nc), probability=node_pb) 
					nc += 1
					node_add += 1
					
		leaf_nodes = [x for x in self.G.nodes_iter() if self.G.out_degree(x) == 0]
		
		for n in leaf_nodes:
			
			## Determine the path between root node and previous time node
			rpath = networkx.all_simple_paths(self.G,'R',n)
			rpath = list(rpath)
			
			if len(rpath) > 1: 
				print('Error in Tree')
			
			if len(rpath) == 0:
				cummulative_decisions = self.G.node['R']['Fixed Variables']
			else:
				all_decs = [self.G.node[nint]['Fixed Variables'] for nint in rpath[0]]
				cummulative_decisions = []
				for itr in all_decs:
					cummulative_decisions += itr 
			
			### Identify what is the next step decisions
			Trial_Constraints = [0 for x in ModelData._data['product'][None]]

			for (i,j,t,d) in cummulative_decisions:
				
				i_idx = ModelData._data['product'][None].index(i)
				j_idx = ModelData._data['trial'][None].index(j)
				
				## Sets the index of the constrained trial for each product. If all are zero then the first trial has not been started
				if j_idx >= Trial_Constraints[i_idx] and d == 1:
					Trial_Constraints[i_idx] = j_idx + 1
					
			### Once we identified the decisions which have been constrained we can look at the last node and determine the next time period in that is differentiable
			
			ct = self.G.node[n]['Active Periods']
			
			### Minimum of trials that can be started
			MTD = 100
			p = 0
			try:
				while p < len(Trial_Constraints):
					pr = ModelData._data['product'][None][p]
					if Trial_Constraints[p] < len(ModelData._data['trial'][None]):
						tr = ModelData._data['trial'][None][Trial_Constraints[p]]
						if ModelData._data['trial_duration'][(pr,tr)] < MTD:
							MTD = ModelData._data['trial_duration'][(pr,tr)]	
					p += 1
			except:
				pdb.set_trace()
			
			### This converts MTD from Time Periods to Time of Realization
			MTD += ct
			
			### If its a pseudonode you need to remove 1 (bc it doesn't really exist)
			if self.G.node[str(n)]['Fixed Variables']==[] and n.startswith('u'):
				MTD = MTD - 1
				
				### we need to add constraints to account for cases where there are no additional pseudo nodes
				self.G.node[str(n)]['Constraints'] = [Trial_Constraints[i] for i in range(len(ModelData._data['product'][None]))]
				
				
			### Active trials in the pipeline
			SR = 100
			
			for (i,j,t,d) in cummulative_decisions:
				if d == 1:
					td = ModelData._data['trial_duration'][(i,j)]
					if t + td > ct:
						if t + td - 1 < SR:
							SR = t + td -1
								
			### Create Pseudo nodes based on the active trials and the trials that can be started
			# The soonest anything that can be realized is the minimum of MTD and SR
			end_ps_node = min(SR,MTD)
			
			### No pseudo nodes!
			if end_ps_node == ct:
				self.G.node[str(n)]['Constraints'] = [Trial_Constraints[i] for i in range(len(ModelData._data['product'][None]))]
			else:
				ap = 0
				while ap < end_ps_node - ct :
					if ap + ct < ModelData._data['time_step'][None][-1]:
						### add pseudo nodes
						self.G.add_node('up' + str(nc))
						self.G.node['up' + str(nc)]['Active Periods'] = ct + 1 + ap
						self.G.node['up' + str(nc)]['Fixed Variables'] = []
						self.G.node['up' + str(nc)]['Constraints'] = [Trial_Constraints[i] for i in range(len(ModelData._data['product'][None]))]
						self.G.node['up' + str(nc)]['Realizations'] = copy.deepcopy(self.G.node[n]['Realizations'])
						if str("T"+ str(ct + ap + 1)) not in self.stage_names:
							self.stage_names.append("T"+ str(ct + 1 + ap))
						if 1 + ap == 1:
							self.G.add_edge(n,'up' + str(nc), probability = 1)
						else:
							self.G.add_edge('up' + str(nc-1),'up' + str(nc), probability = 1)
						nc += 1
					ap += 1
						
	def add_scenarios(self, ModelData):		
		### Generate scenarios
		self.stage_names.append("Ts")
		
		### Define number of products and number of trials
		NP = len(ModelData._data['product'][None])
			
		NT = len(ModelData._data['trial'][None])
		
		### Create Scenario Nodes
		Outcomes = itertools.product(range(NT + 1), repeat = NP)
		Outcomes = list(Outcomes)
	
		### Generate Leaf Nodes
		leaf_nodes = [x for x in self.G.nodes_iter() if self.G.out_degree(x) == 0]
		
		tp = {}
		pn = 0
		for s in Outcomes:
			
			### Determine which leaf node to attach to
			matched = False	
			cntr = 0	
			while not matched and cntr < len(leaf_nodes):
				for n in leaf_nodes:
					match = True
					for (i,j,r) in self.G.node[n]['Realizations']:
							
						### Calculate whether realizations match
						if r == 1:
							
							### The trial (i,j) passed, the outcome must be greater than or equal to 
							idx = ModelData._data['product'][None].index(i)
								
							if s[idx] >= j:
								### This is a match
								pass
							else:
								match = False
						else:
							### The trial (i,j) passed, the outcome must be greater than or equal to 
							idx = ModelData._data['product'][None].index(i)
								
							if s[idx] < j:
								### This is a match
								pass
							else: 
								match = False
									
						### made it through all the realizations
					if match:	
						match_node = n
						matched = True
					cntr += 1
					
			if not matched:
				pdb.set_trace()
				
			### Calculate Conditional Probability (the absolute is equivalent to the conditional for scenarios) 
			p = 1
			i = 0
			while i < len(s):
				if s[i] == len(ModelData._data['trial'][None]):
					for j in ModelData._data['trial'][None]:
						coords = (ModelData._data['product'][None][i],j)
						p *= ModelData._data['probability'][coords]
				else:
					for j in ModelData._data['trial'][None]:
						tindex = ModelData._data['trial'][None].index(j)
						coords = (ModelData._data['product'][None][i],ModelData._data['trial'][None][tindex])
						if tindex < s[i]:
							p *= ModelData._data['probability'][coords]
						elif tindex == s[i]:
							p *= (1 - ModelData._data['probability'][coords])
				i += 1
			
			if match_node in tp:
				tp[match_node] += p
			else:
				tp[match_node] = p
			
			### Check to make sure that there are the appropriate number of stages, if there are not add ps nodes
			
			## Determine the path between root node and leaf node
			rpath = networkx.all_simple_paths(self.G,'R',match_node)
			rpath = list(rpath)
			
			if match_node != 'R':
				if len(rpath) == 0:
					### Well no path exists
					pdb.set_trace()
				
			pn_added = False
			pn = self.G.node[match_node]['Active Periods'] + 1
			if match_node != 'R':
				if len(rpath[0]) + 1 < len(self.stage_names):
			
					### There need to be ps stages added bc attaching scenarios will cause solver failure
					nodes_needed = len(self.stage_names) - len(rpath[0]) - 1
					pn_added = True
					for nad in range(nodes_needed):
						self.G.add_node('pn_'+ str(pn)+ '_' + str(Outcomes.index(s)))
						if nad == 0:
							self.G.add_edge(match_node,'pn_'+ str(pn)+ '_' + str(Outcomes.index(s)),probability=p)
							self.G.node['pn_'+ str(pn)+ '_' + str(Outcomes.index(s))]['Active Periods']= self.G.node[match_node]['Active Periods'] + 1
							pn += 1
						else:
							self.G.add_edge('pn_'+ str(pn-1)+ '_' + str(Outcomes.index(s)),'pn_'+ str(pn)+ '_' + str(Outcomes.index(s)),probability=1)
							self.G.node['pn_'+ str(pn)+ '_' + str(Outcomes.index(s))]['Active Periods']= self.G.node['pn_'+ str(pn-1)+ '_' + str(Outcomes.index(s))]['Active Periods'] + 1
							pn += 1
					
			### Attach Node
			if pn_added == False:
				self.G.add_node('s' + str(Outcomes.index(s)))
				self.G.node['s' + str(Outcomes.index(s))]['Outcome'] = list(s)
				self.G.add_edge(match_node,'s' + str(Outcomes.index(s)),probability=p)
				self.G.node['s' + str(Outcomes.index(s))]['Active Periods']= self.G.node[match_node]['Active Periods'] + 1
				self.G.node['s' + str(Outcomes.index(s))]['PB'] = p
			else:
				self.G.add_node('s' + str(Outcomes.index(s)))
				self.G.node['s' + str(Outcomes.index(s))]['Outcome'] = list(s)
				self.G.node['s' + str(Outcomes.index(s))]['Active Periods']= self.G.node['pn_'+ str(pn-1)+ '_' + str(Outcomes.index(s))]['Active Periods'] + 1
				self.G.node['s' + str(Outcomes.index(s))]['PB'] = p
				self.G.add_edge('pn_'+ str(pn-1)+ '_' + str(Outcomes.index(s)),'s' + str(Outcomes.index(s)),probability=1)
			
				
		### Update to conditional probabilities 
		for n in leaf_nodes:
			for node_scenario in self.G.edge[n]:
				self.G.edge[n][node_scenario]['probability'] = self.G.edge[n][node_scenario]['probability']/ tp[n]
			

