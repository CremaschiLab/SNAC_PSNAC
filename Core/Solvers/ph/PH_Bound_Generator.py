import sys
import os
from pyomo.environ import *
from pyomo.pysp.scenariotree.tree_structure_model import ScenarioTreeModelFromNetworkX
from pyomo.pysp.phinit import exec_runph as exec_runph
from pyomo.pysp.phinit import main as ph_instance
from Core.Solvers.PH.Decision_Tree_Object import Multistage_PYSP_Object as DT
import networkx
import pdb
from pyomo.pysp.scenariotree.manager_solver import \
    ScenarioTreeManagerSolverClientSerial
import builtins
import Core.Solvers.MTSSP.M2S_item as M2S_item
import copy
from pyomo.pysp.util.misc import (launch_command,
                                  load_extensions)
import csv

	
def Bound_Generator(MD,dec, output_directory):
	
	# declare the number of scenarios over which to construct a simple
	# two-stage scenario tree
	
	set_global(MD,dec)
	
	## Determine the Path to this file
	thisfile = os.path.abspath(__file__)
	
	manager_type = ScenarioTreeManagerSolverClientSerial
	options = ScenarioTreeManagerSolverClientSerial.register_options()
	
	### Detect and use the pysp_instance_creation_callback() function
	options.model_location = thisfile
	
	### Determine the location of the scenario tree file, None indicates 
	## that the scenario tree is defined in the modelfile
	options.scenario_tree_location = None
	
	dir_loc = str(thisfile)
	scen_tree_loc = '/home/che_h2/bzc0043/Dropbox/Katie_Research_Team_Folder/Code/SPHeuristic_1.6.5/Core/Solvers/PH/Scenario_Tree_Structure.py'
	
	save_file = 'ph_output'
	output_file = os.path.join(output_directory,save_file)
	
	
	results = ph_instance(('-m',dir_loc, '-i',scen_tree_loc,'--default-rho','1','--solution-writer=pyomo.pysp.plugins.csvsolutionwriter','--traceback', '--max-iterations=200'))
		
	# Fetch bound info from csv	(this should not have to be done!!!!!!)
	csvfile = open('ph_StageCostDetail.csv') 
		
	ph_reader = csv.reader(csvfile)
	
	obj_val = 0
	
	for row in ph_reader:
		### Get probability of scenario row[2]
		sname = row[2].replace(" ","")
		innode = builtins.Tree_Graph_GLOBAL.G.in_edges(sname)
		innode = innode[0][0]
		pb = builtins.Tree_Graph_GLOBAL.G.node[sname]['PB']
		
		### Multiply probability by row[5] and add to obj_val
		obj_val += pb * float(row[5])
		
	csvfile.close()
	
	obj_val = -1 * obj_val
	
	return obj_val

def return_fixed_values(i,j,t,nd):
	for (ip,jp,tp,decp) in builtins.Tree_Graph_GLOBAL.G.node[nd]['Fixed Variables']:
		if i==ip and j==jp and t == tp:
			return decp 
	print('Error in locating fixed variable')
	exit()
	
def set_global(MD, dec):
	builtins.MD_GLOBAL = MD
	builtins.dec_GLOBAL = dec
	return

def pysp_instance_creation_callback(scenario_name, node_names):
	### Build the Scenario Model
	m = ConcreteModel()
		
	##Sets##
	m.I = Set(initialize=builtins.MD_GLOBAL._data['product'][None])
	m.J = Set(initialize=builtins.MD_GLOBAL._data['trial'][None])
	m.T = Set(initialize=builtins.MD_GLOBAL._data['time_step'][None])
	m.R = Set(initialize=builtins.MD_GLOBAL._data['resource_type'][None])
	
		### Fixed Parameters###
	m.Resource_Max = {}
	for r in builtins.MD_GLOBAL._data['resource_type'][None]:
		m.Resource_Max[r] = builtins.MD_GLOBAL._data['max_resource'][(r,)]
	
	## Set Discount Values
	m.GammaL = {}
	for items in builtins.MD_GLOBAL._data['gammaL']:
		m.GammaL[items[0]] = builtins.MD_GLOBAL._data['gammaL'][items]
	
	m.GammaD = {}	
	for items in builtins.MD_GLOBAL._data['gammaD']:
		m.GammaD[items[0]] = builtins.MD_GLOBAL._data['gammaD'][items]
	
	## Set Maximum Revenue
	m.Revenue_Max = {}	
	for items in builtins.MD_GLOBAL._data['maximum_revenue']:
		m.Revenue_Max[items[0]] = builtins.MD_GLOBAL._data['maximum_revenue'][items]
	
	m.Duration = builtins.MD_GLOBAL._data['trial_duration'] 
	m.Trial_Cost = builtins.MD_GLOBAL._data['trial_cost'] 
	m.Resources_Required = builtins.MD_GLOBAL._data['resource_requirement']
	
	
	m.LastTimeStep = builtins.MD_GLOBAL._data['time_step'][None][-1]
	m.Last_Trial = builtins.MD_GLOBAL._data['trial'][None][-1]
	
	## Calculate running rev
	m.RR = M2S_item.calc_rr(m.Revenue_Max,m.GammaL,m.Duration, m.I, m.J, m.T)
		
	##Calculate open rev  
	m.OR = M2S_item.calc_openrev(m.Revenue_Max,m.GammaL,m.Duration, m.I, m.J, m.T, m.LastTimeStep)
	
	##Calculate Discounting Factor
	m.DF = M2S_item.calc_discounting_factor(m.Revenue_Max,m.GammaL,m.Trial_Cost, m.I, m.J, m.LastTimeStep)
	
	### Scenario Specific Parameters###
	m.Success = {}
	
	for i in builtins.MD_GLOBAL._data['product'][None]:
		if builtins.Tree_Graph_GLOBAL.G.node[scenario_name]['Outcome'][builtins.MD_GLOBAL._data['product'][None].index(i)] == len(m.J):
			m.Success[i] = 1
		else:
			m.Success[i] = 0	
	
	def FS_init(node):
		rval = []
		t = builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']
		tc = 0
		
		### If a node isn't a isn't a 'up' node then all variables are first stage
		if node.startswith('up'):
			
			## Find Last Pseudonode
			for nds in reversed(node_names):
				if 'Constraints' in builtins.Tree_Graph_GLOBAL.G.node[nds]:
					last_pseudo = nds
			try:	
				cvars = builtins.Tree_Graph_GLOBAL.G.node[last_pseudo]['Constraints']
			except:
				### Error thrown relates to not having a last pseudo
				pdb.set_trace()
				
			while tc < len(cvars):
				for trial in builtins.MD_GLOBAL._data['trial'][None]:
					if builtins.MD_GLOBAL._data['trial'][None].index(trial) <= cvars[tc]:
						if (builtins.MD_GLOBAL._data['product'][None][tc],trial,t,0) not in builtins.Tree_Graph_GLOBAL.G.node[node]['Fixed Variables'] and (builtins.MD_GLOBAL._data['product'][None][tc],trial,t,1) not in builtins.Tree_Graph_GLOBAL.G.node[node]['Fixed Variables']:
							rval.append((builtins.MD_GLOBAL._data['product'][None][tc],trial,t))
			
				tc += 1
			return tuple(rval)
		else:
			for i in builtins.MD_GLOBAL._data['product'][None]:
				for j in builtins.MD_GLOBAL._data['trial'][None]:
					if (i,j,t,0) not in builtins.Tree_Graph_GLOBAL.G.node[node]['Fixed Variables'] and (i,j,t,1) not in builtins.Tree_Graph_GLOBAL.G.node[node]['Fixed Variables']:
						rval.append((i,j,t))
			return tuple(rval)
		
	def Fixed_init(node):
		rval = []
		for (i,j,t,r) in builtins.Tree_Graph_GLOBAL.G.node[node]['Fixed Variables']:
			rval.append((i,j,t))
		return tuple(rval)	
		
	def SS_init(node):
		rval = []
		t = builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']
		tc = 0
		### If a node isn't a isn't a 'up' node then all variables are first stage
		if node.startswith('up'):	
			## Find Last Pseudonode
			for nds in reversed(node_names):
				if 'Constraints' in builtins.Tree_Graph_GLOBAL.G.node[nds]:
					last_pseudo = nds
					
			cvars = builtins.Tree_Graph_GLOBAL.G.node[last_pseudo]['Constraints']
				
			while tc < len(cvars):
				if cvars[tc] < len(builtins.MD_GLOBAL._data['trial'][None]):
					trl = builtins.MD_GLOBAL._data['trial'][None][cvars[tc]]
					for trial in builtins.MD_GLOBAL._data['trial'][None]:
						if builtins.MD_GLOBAL._data['trial'][None].index(trial) > cvars[tc]:
							if (builtins.MD_GLOBAL._data['product'][None][tc],trial,t,0) not in builtins.Tree_Graph_GLOBAL.G.node[node]['Fixed Variables'] and (builtins.MD_GLOBAL._data['product'][None][tc],trial,t,1) not in builtins.Tree_Graph_GLOBAL.G.node[node]['Fixed Variables']:
								rval.append((builtins.MD_GLOBAL._data['product'][None][tc],trial,t))
				tc += 1
			return tuple(rval)
		else:
			return(tuple(rval))
	
	def scen_init(node):
		### all decisions after t
		tend = builtins.Tree_Graph_GLOBAL.G.node[node_names[-2]]['Active Periods']
		rval = []
		for t in builtins.MD_GLOBAL._data['time_step'][None]:
			if t > tend:
				for i in builtins.MD_GLOBAL._data['product'][None]:
					for j in builtins.MD_GLOBAL._data['trial'][None]:
						rval.append((i,j,t))
		return tuple(rval)
	
	def pseudo_scen_init(node):
		rval = []
		for i in builtins.MD_GLOBAL._data['product'][None]:
			for j in builtins.MD_GLOBAL._data['trial'][None]:
				rval.append((i,j,builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']))
		return rval
				
	nodeallset = {}
	nodevar = {}	
	nodeset = {}
	nodeallvar = {}
	for nn in node_names:
		nodeset[nn] = {}
		nodevar[nn] = {}
		if nn.startswith('s'):
			m.add_component(str(nn) + 'SetFixed', Set(within = m.I * m.J * m.T, initialize = set()))
			m.add_component(str(nn) + 'SetX', Set(within = m.I * m.J * m.T, initialize = set()))
			m.add_component(str(nn) + 'SetXX', Set(within = m.I * m.J * m.T, initialize = scen_init(nn)))
		elif nn.startswith('pn'):
			m.add_component(str(nn) + 'SetFixed', Set(within = m.I * m.J * m.T, initialize = set()))
			m.add_component(str(nn) + 'SetX', Set(within = m.I * m.J * m.T, initialize = set()))
			m.add_component(str(nn) + 'SetXX', Set(within = m.I * m.J * m.T, initialize = pseudo_scen_init(nn)))
		else:
			m.add_component(str(nn) + 'SetFixed', Set(within = m.I * m.J * m.T, initialize = Fixed_init(nn)))
			m.add_component(str(nn) + 'SetX', Set(within = m.I * m.J * m.T, initialize = FS_init(nn)))
			m.add_component(str(nn) + 'SetXX', Set(within = m.I * m.J * m.T, initialize = SS_init(nn)))
	
	
		
	for nn in node_names:
		for comp in list(m.component_objects()):
			if comp._name.startswith(str(nn) + 'Set'):
				if not comp._name.endswith('_domain') and not comp._name.endswith('_domain_index_0'):
					nodeset[nn][str(comp._name)]= comp
					if nn.startswith('s') or nn.startswith('pn'):
						if comp._name.startswith(str(nn) + 'SetXX'):
							m.add_component(str(nn) + 'VarXX', Var(comp,bounds=(0,1)))
					else:
						if comp._name.startswith(str(nn) + 'SetXX'):
							m.add_component(str(nn) + 'VarXX', Var(comp,bounds=(0,1)))
						elif comp._name.startswith(str(nn) + 'SetX'):
							m.add_component(str(nn) + 'VarX', Var(comp,bounds=(0,1)))
	for nn in nodeset:
		setAll = set()
		nodeallset[nn] = {}
		for sets in nodeset[nn]:
			setNew = setAll.union(set(nodeset[nn][sets].value))
			setAll = setNew
		m.add_component(str(nn) + '_setALL', Set(within=m.I*m.J*m.T,initialize = copy.deepcopy(setAll)))
					
	for nn in node_names:
		for comp in list(m.component_objects()):		
			if comp._name == str(nn) + '_setALL': 				
				nodeallset[nn][str(comp._name)] = comp
				m.add_component(str(nn)+'VarY',Var(comp,bounds=(0,1)))
				m.add_component(str(nn) + 'VarZ', Var(comp,bounds=(0,1)))
	
	for nn in node_names:
		nodeallvar[nn] = {}
		for comp in list(m.component_objects()):
			
			if comp._name == str(nn) + 'VarY':
				nodeallvar[nn]['Y'] = comp
			elif comp._name == str(nn) + 'VarZ':
				nodeallvar[nn]['Z'] = comp
								
	for nn in node_names:
		for comp in list(m.component_objects()):
			if comp._name.startswith(str(nn) + 'VarX') or comp._name.startswith(str(nn) + 'VarXX') :
				nodevar[nn][str(comp._name)] = comp
	
	if builtins.Tree_Graph_GLOBAL.debug == True:
		pdb.set_trace()
	
	def construct_stage_cost(model,ndn):
	
		cst = 0
		### Cost to start trials for Fixed Variables, Constrained Variables, and Second Stage Variables
		#pdb.set_trace()
		for (i,j,t) in nodeallset[ndn][str(ndn)+'_setALL']:
			
			if (i,j,t) in nodeset[ndn][str(ndn)+'SetX']:
				### Cost to start trials
				cst += -(1 - (0.025 * (t - 1))) * model.Trial_Cost[i,j] * nodevar[ndn][str(ndn)+'VarX'][i,j,t]
				
				### Revenue From Starting Final Trials
				if j == m.Last_Trial:
					cst += m.Success[i]*(m.Revenue_Max[i] * nodevar[ndn][str(ndn)+'VarX'][i,j,t])
					
				
				### Losses for the time to develop drug
				if j == m.Last_Trial:
					cst -= m.Success[i]*m.GammaL[i]*(t + m.Duration[i, j]) * nodevar[ndn][str(ndn)+'VarX'][i,j,t]
					
				
				### Revenue not yet realized
				if j < m.Last_Trial and t > m.LastTimeStep - m.Duration[i,j]:
					cst += m.Success[i] * m.RR[i,j,t] * m.DF[i,j+1]*nodevar[ndn][str(ndn)+'VarX'][i,j,t]
					
		
			elif (i,j,t) in nodeset[ndn][str(ndn)+'SetFixed']:
				
				### Cost to start trials
				cst += -(1 - (0.025 * (t - 1))) * model.Trial_Cost[i,j] * return_fixed_values(i,j,t,ndn) 
				
				### Revenue From Starting Final Trials
				if j == m.Last_Trial:
					cst += m.Success[i]*(m.Revenue_Max[i] * return_fixed_values(i,j,t,ndn)) 
				
				### Losses for the time to develop drug
				if j == m.Last_Trial:
					cst += -m.Success[i]*m.GammaL[i]*(t + m.Duration[i, j]) * return_fixed_values(i,j,t,ndn)
					
					
				
				### Revenue not yet realized
				if j < m.Last_Trial and t > m.LastTimeStep - m.Duration[i,j]:
					cst += m.Success[i] * m.RR[i,j,t] * m.DF[i,j+1] * return_fixed_values(i,j,t,ndn)
					 
			else:
				
				### Cost to start trials
				cst += -(1 - (0.025 * (t - 1))) * model.Trial_Cost[i,j] * nodevar[ndn][str(ndn)+'VarXX'][i,j,t]
				
				### Revenue From Starting Final Trials
				if j == m.Last_Trial:
					cst += m.Success[i]*(m.Revenue_Max[i] * nodevar[ndn][str(ndn)+'VarXX'][i,j,t])
				
				
				### Losses for the time to develop drug
				if j == m.Last_Trial:
					cst -= m.Success[i]*m.GammaL[i]*(t + m.Duration[i, j]) * nodevar[ndn][str(ndn)+'VarXX'][i,j,t] 
					
				
				### Revenue not yet realized
				if j < m.Last_Trial and t > m.LastTimeStep - m.Duration[i,j]:
					cst += m.Success[i] * m.RR[i,j,t] * m.DF[i,j+1]*nodevar[ndn][str(ndn)+'VarXX'][i,j,t]
		
		if ndn.startswith('s'):
			tstart = builtins.Tree_Graph_GLOBAL.G.node[ndn]['Active Periods']
			
			### Losses For Decisions Not Made
			cst += sum(-m.Success[i]*m.GammaD[i] * sum(nodeallvar[ndn]['Z'][i,j,t] for j in m.J if j > 1) for i in m.I for t in m.T if t >= tstart)
			
			if m.LastTimeStep >= tstart:
				cst += sum(m.Success[i] * m.OR[i,j] * m.DF[i,j]* nodeallvar[ndn]['Z'][i,j,m.LastTimeStep] for i in m.I for j in m.J)
				
		else:
			t = builtins.Tree_Graph_GLOBAL.G.node[ndn]['Active Periods']	
		
			### Losses For Decisions Not Made
			cst += sum(-m.Success[i]*m.GammaD[i] * sum(nodeallvar[ndn]['Z'][i,j,t] for j in m.J if j > 1) for i in m.I)
		
			
			### Revenue Yet to Be Realized
			if t == m.LastTimeStep:
				cst += sum(m.Success[i] * m.OR[i,j] * m.DF[i,j]* nodeallvar[ndn]['Z'][i,j,t] for i in m.I for j in m.J)
				
		return -cst
	
	## Add stage cost for each node
	stage_covered = []
	for node in node_names:
		if node.startswith('s'):
			### This is the scenario stage
			match_node = scenario_name
			stage_covered.append('Ts')
			m.add_component('TsStageCost', Expression(expr=construct_stage_cost(m, match_node)))
		
		else:	
			### Determine stage
			tstage = builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']
			
			stage_name = 'T' + str(tstage)
			stage_covered.append(stage_name)
			
			try:
				m.add_component(str(stage_name) + 'StageCost', Expression(expr=construct_stage_cost(m, node)))	
			except Exception as e:
				print(e)
				pdb.set_trace()
				m.add_component(str(stage_name) + 'StageCost', Expression(expr=construct_stage_cost(m, node)))
					
	for stage in builtins.Tree_Graph_GLOBAL.stage_names:
		if stage not in stage_covered:
			m.add_component(str(stage) + 'StageCost', Expression(expr=0.0))
				

	m.obj = Objective(expr = sum(comp for comp in list(m.component_objects()) if comp._name.endswith('StageCost')))
	
	### Constraint-- Ensures resources are managed correctly
	def Resource_Constraint_rule(m,r,t):
				
		LHS = 0
		for i in m.I:
			for j in m.J:
				for tprime in m.T:
					
					if tprime > (t - m.Duration[i,j]) and tprime <= t:
						##identify node
						match_node = [n for n in node_names if builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods'] == tprime]
						
						if len(match_node) > 1 or len(match_node) == 0:
							if len(match_node) == 0:
								if tprime > builtins.Tree_Graph_GLOBAL.G.node[node_names[-1]]['Active Periods']:
									match_node = node_names[-1]
								else:
									print('You done bunked it')
						else:
							match_node = match_node[0]
		
						if match_node.startswith('s') or match_node.startswith('pn') :
							LHS += m.Resources_Required[i,j,r] * nodevar[match_node][str(match_node)+'VarXX'][i,j,tprime]
						else:
							if (i,j,tprime) in nodeset[match_node][str(match_node)+'SetX']:
								LHS += m.Resources_Required[i,j,r] * nodevar[match_node][str(match_node)+'VarX'][i,j,tprime]
							
							elif (i,j,tprime) in nodeset[match_node][str(match_node)+'SetFixed']:
								LHS += m.Resources_Required[i,j,r] * return_fixed_values(i,j,tprime,match_node)  
								
							else:
								LHS += m.Resources_Required[i,j,r] * nodevar[match_node][str(match_node)+'VarXX'][i,j,tprime]
						
		return LHS <= m.Resource_Max[r]	
			
	m.Resource_Constraint = Constraint(m.R, m.T, rule= Resource_Constraint_rule)

		
	### Constraint-- Must Start Previous Trial Before it can be finished
	def Constraint_4_rule(m,i,j,t):
		if j > 1:
			previous_trial = j-1
			LHS = 0
			for tprime in m.T:
				
				if tprime <= t:
					### Determine the node
					match_node = [n for n in node_names if builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods'] == tprime]
						
					if len(match_node) > 1 or len(match_node) == 0:
						if len(match_node) == 0:
							if tprime >= builtins.Tree_Graph_GLOBAL.G.node[node_names[-1]]['Active Periods']:
								match_node = node_names[-1]
							else:
								print('You done bunked it')
					else:
						match_node = match_node[0]
			
					if match_node.startswith('s') or match_node.startswith('pn'):
						LHS += nodevar[match_node][str(match_node)+'VarXX'][i,j,tprime]
					else:
						if (i,j,tprime) in nodeset[match_node][str(match_node)+'SetX']:
							LHS += nodevar[match_node][str(match_node)+'VarX'][i,j,tprime]
							
						elif (i,j,tprime) in nodeset[match_node][str(match_node)+'SetFixed']:
							LHS += return_fixed_values(i,j,tprime,match_node)  
								
						else:
							LHS += nodevar[match_node][str(match_node)+'VarXX'][i,j,tprime]
			### Determine node for t
			match_node = [n for n in node_names if builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods'] == t]
						
			if len(match_node) > 1 or len(match_node) == 0:
				if len(match_node) == 0:
					if tprime >= builtins.Tree_Graph_GLOBAL.G.node[node_names[-1]]['Active Periods']:
						match_node = node_names[-1]
					else:
						print('You done bunked it')
			else:
				match_node = match_node[0]
				
			decY = nodeallvar[match_node]['Y'][i,j-1,t]
			return LHS <= decY
		else:
			return Constraint.Skip
			
	m.Constraint_4 = Constraint(m.I, m.J, m.T, rule=Constraint_4_rule)
	
	### Constraint--You can only start each trial once
	def Constraint_3_rule(model,i,j):
		LHS = 0
		for t in m.T:
			### Determine the node
			match_node = [n for n in node_names if builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods'] == t]
						
			if len(match_node) > 1 or len(match_node) == 0:
				if len(match_node) == 0:
					if t >= builtins.Tree_Graph_GLOBAL.G.node[node_names[-1]]['Active Periods']:
						match_node = node_names[-1]
					else:
						print('You done bunked it')
			else:
				match_node = match_node[0]
							
			if match_node.startswith('s') or match_node.startswith('pn'):
				LHS += nodevar[match_node][str(match_node)+'VarXX'][i,j,t]
			else:
				if (i,j,t) in nodeset[match_node][str(match_node)+'SetX']:
					LHS += nodevar[match_node][str(match_node)+'VarX'][i,j,t]
							
				elif (i,j,t) in nodeset[match_node][str(match_node)+'SetFixed']:
					LHS += return_fixed_values(i,j,t,match_node) 
								
				else:
					LHS += nodevar[match_node][str(match_node)+'VarXX'][i,j,t]
		return LHS <= 1
	m.Constraint_3 = Constraint(m.I,m.J,rule=Constraint_3_rule)
	
	def Calculate_Z_rule(m,i,j,t):
		
		### Determine the node
		match_node = [n for n in node_names if builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods'] == t]
						
		if len(match_node) > 1 or len(match_node) == 0:
			if len(match_node) == 0:
				if t >= builtins.Tree_Graph_GLOBAL.G.node[node_names[-1]]['Active Periods']:
					match_node = node_names[-1]
				else:
					print('You done bunked it')
		else:
			match_node = match_node[0]
		
		### Determine X,Z,Zpd and Xpd	
		if match_node.startswith('s') or match_node.startswith('pn'):
			decX = nodevar[match_node][str(match_node)+'VarXX'][i,j,t]
			decZ = nodeallvar[match_node]['Z'][i,j,t]
			
			### Find previous nodes for Xpd
			if j-1 > 0:
				if t-m.Duration[i,j-1] > 0:
					if t-m.Duration[i,j-1] >= builtins.Tree_Graph_GLOBAL.G.node[match_node]['Active Periods']:
						
						### This means that Xpd is inside the set for S
						decXpd = nodevar[match_node][str(match_node)+'VarXX'][i,j-1,t-m.Duration[i,j-1]]
					else:
						for node in node_names:
							if t-m.Duration[i,j-1] == builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']:
								### The node is the node now check FD, FSD or SSD
								if (i,j-1,t-m.Duration[i,j-1]) in nodeset[node][str(node)+'SetX']:
									decXpd = nodevar[node][str(node)+'VarX'][i,j-1,t-m.Duration[i,j-1]]
								elif (i,j-1,t-m.Duration[i,j-1]) in nodeset[node][str(node)+'SetFixed']:
									decXpd = return_fixed_values(i,j-1,t-m.Duration[i,j-1],node) 
								else:
									decXpd = nodevar[node][str(node)+'VarXX'][i,j-1,t-m.Duration[i,j-1]]
			if t-1 > 0:
				if t-1 >= builtins.Tree_Graph_GLOBAL.G.node[match_node]['Active Periods']:
					
					### This means that Zpt is inside the set for S
					decZpt = nodeallvar[match_node]['Z'][i,j,t-1]
				else:
					for node in node_names:
						if t-1 == builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']:
							
							### The node is Z
							decZpt = nodeallvar[node]['Z'][i,j,t-1]
					
		else:
			if (i,j,t) in nodeset[match_node][str(match_node)+'SetX']:
				decX = nodevar[match_node][str(match_node)+'VarX'][i,j,t]
			elif (i,j,t) in nodeset[match_node][str(match_node)+'SetFixed']:
				decX = return_fixed_values(i,j,t,match_node)
			else:
				decX = nodevar[match_node][str(match_node)+'VarXX'][i,j,t]
			
			decZ = nodeallvar[match_node]['Z'][i,j,t]
			
			
			for node in node_names:
				if j-1 > 0:	
					if t-m.Duration[i,j-1] > 0:
						if t-m.Duration[i,j-1] == builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']:
							### The node is the node now check FD, FSD or SSD
							if (i,j-1,t-m.Duration[i,j-1]) in nodeset[node][str(node)+'SetX']:
								decXpd = nodevar[node][str(node)+'VarX'][i,j-1,t-m.Duration[i,j-1]]
							elif (i,j-1,t-m.Duration[i,j-1]) in nodeset[node][str(node)+'SetFixed']:
								decXpd = return_fixed_values(i,j-1,t-m.Duration[i,j-1],node) 
							else:
								decXpd = nodevar[node][str(node)+'VarXX'][i,j-1,t-m.Duration[i,j-1]]
				if t - 1 > 0:
					if t-1 == builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']:
						### The node is Z
						decZpt = nodeallvar[node]['Z'][i,j,t-1]		
		
		decZ1 = nodeallvar['R']['Z'][i,j,1]		
		if j == 1:
			if t > 1:
				return decZ == decZpt - decX
			else:
				return decZ1 == 1 - decX
		else:
			if t- m.Duration[i,j-1] >0 and t>1:
				return decZ == decZpt + decXpd - decX
			elif t-m.Duration[i,j-1] <= 0 and t > 1:
					return decZ == decZpt - decX
			elif t==1:
				return decZ == -decX
	
	m.Calculate_Z = Constraint(m.I,m.J,m.T, rule= Calculate_Z_rule)

	
	### Constraint-- The trial is complete if the trial was started a duration ago
	def Trial_Finish_rule(m,i,j,t):
		past_duration = t - m.Duration[i,j]
		
		### Determine the node
		match_node = [n for n in node_names if builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods'] == t]
						
		if len(match_node) > 1 or len(match_node) == 0:
			if len(match_node) == 0:
				if t >= builtins.Tree_Graph_GLOBAL.G.node[node_names[-1]]['Active Periods']:
					match_node = node_names[-1]
				else:
					print('You done bunked it')
		else:
			match_node = match_node[0]
		
		##decXpd and decYpt
		### Find previous nodes for Xpd
			
		if t-m.Duration[i,j] > 0:
			if t-m.Duration[i,j] >= builtins.Tree_Graph_GLOBAL.G.node[match_node]['Active Periods']:
					
					### This means that Xpd is inside the set for S
				decXpd = nodevar[match_node][str(match_node)+'VarXX'][i,j,t-m.Duration[i,j]]
			else:
				for node in node_names:
					if t-m.Duration[i,j] == builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']:
						### The node is the node now check FD, FSD or SSD
						if node.startswith('pn'):
							decXpd = nodevar[node][str(node)+'VarXX'][i,j,t-m.Duration[i,j]]
						else:
							if (i,j,t-m.Duration[i,j]) in nodeset[node][str(node)+'SetX']:
								decXpd = nodevar[node][str(node)+'VarX'][i,j,t-m.Duration[i,j]]
							elif (i,j,t-m.Duration[i,j]) in nodeset[node][str(node)+'SetFixed']:
								decXpd = return_fixed_values(i,j,t-m.Duration[i,j],node) 
							else:
								decXpd = nodevar[node][str(node)+'VarXX'][i,j,t-m.Duration[i,j]]
		if t-1 > 0:
			if t-1 >= builtins.Tree_Graph_GLOBAL.G.node[match_node]['Active Periods']:
				
				### This means that Zpt is inside the set for S
				decYpt = nodeallvar[match_node]['Y'][i,j,t-1]
			else:
				for node in node_names:
					if t-1 == builtins.Tree_Graph_GLOBAL.G.node[node]['Active Periods']:
						
						### The node is Z
						decYpt = nodeallvar[node]['Y'][i,j,t-1]
		
		if t==1:
			return nodeallvar['R']['Y'][i,j,1] == 0
		elif  t > 1 and  past_duration < 1 :
			return nodeallvar[match_node]['Y'][i,j,t] == decYpt
		else:		
			return nodeallvar[match_node]['Y'][i,j,t] == decYpt + decXpd 
	
			
	m.Trial_Finish = Constraint(m.I, m.J, m.T, rule=Trial_Finish_rule)
	
	def trial_start_rule(m,i,j,t):
		jp = builtins.MD_GLOBAL._data['trial'][None].index(j)
		earliest_start = sum(m.Duration[i,jt] for jt in builtins.MD_GLOBAL._data['trial'][None] if builtins.MD_GLOBAL._data['trial'][None].index(jt) < jp)
		
		### Determine the node
		match_node = [n for n in node_names if builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods'] == t]
						
		if len(match_node) > 1 or len(match_node) == 0:
			if len(match_node) == 0:
				if t >= builtins.Tree_Graph_GLOBAL.G.node[node_names[-1]]['Active Periods']:
					match_node = node_names[-1]
				else:
					print('You done bunked it')
		else:
			match_node = match_node[0]
		if t < earliest_start:	
			if match_node.startswith('s') or match_node.startswith('pn'):
				return nodevar[match_node][str(match_node)+'VarXX'][i,j,t] == 0
			else:
				if (i,j,t) in nodeset[match_node][str(match_node)+'SetX']:
					return nodevar[match_node][str(match_node)+'VarX'][i,j,t] == 0
				elif (i,j,t) in nodeset[match_node][str(match_node)+'SetFixed']:
					return Constraint.Skip
				else:
					return nodevar[match_node][str(match_node)+'VarXX'][i,j,t] == 0
		else:
			return Constraint.Skip
	m.Trial_Start = Constraint(m.I,m.J,m.T, rule=trial_start_rule)
	
	return m

	
