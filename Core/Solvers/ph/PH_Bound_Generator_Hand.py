import sys
import os
import math
import multiprocessing as mp
import itertools
from pyomo.environ import *
from Core.Solvers.PH.Decision_Tree_Object import Multistage_PYSP_Object as DT
import networkx
import pdb
import builtins
import Core.Solvers.MTSSP.M2S_item as M2S_item
import copy


	
def Bound_Generator(MD,dec, output_directory):
	
	# declare the number of scenarios over which to construct a simple
	# two-stage scenario tree
	
	set_global(MD,dec)

	## build a decision tree object
	Tree_Graph = DT(dec_GLOBAL,MD_GLOBAL)
	builtins.Tree_Graph_GLOBAL = Tree_Graph
	
	### Define Constants
	rho = 1
	itr = 0
	alpha = 100
	tol = .00001
	max_iterations = 300
	k=0
	
	### For all scenarios we need to 
	scenarios = [x for x in builtins.Tree_Graph_GLOBAL.G.nodes_iter() if builtins.Tree_Graph_GLOBAL.G.out_degree(x) == 0]
	
	### Set MP options
	process_count = mp.cpu_count()
	ng = math.floor(mp.cpu_count()/3)
	np = int(math.ceil(len(scenarios) / (ng)))
	setlist = _grouper(scenarios,np)
	
	### Create Pool
	pool = mp.Pool(ng)
	
	if dec == {('Drug2', 1, 0, 1): [set()], ('Drug3', 1, 0, 0): [set()], ('Drug3', 1, 1, 1): [set()], ('Drug1', 1, 0, 0): [set()], ('Drug2', 2, 2, 1): [{('Drug2', 1, 1)}]}:
		pdb.set_trace()
		for ss in setlist:
			r1, r2 = fs_PH_solve(ss,rho)
	
	results = [pool.apply_async(fs_PH_solve, args=(ss,rho)) for ss in setlist]
	pool.close()
	pool.join()
	
	
	### Merge Results
	
	### Xbar needs to be combine
	xbar = {}
	for n in builtins.Tree_Graph_GLOBAL.G.node:
		if n.startswith('s'):
			pass
		else:
			if n not in xbar:
				xbar[n] = {}
			for r in results:
				try:
					if n in r._value[1]:
						for (i,j,t) in r._value[1][n]:
							if (i,j,t) in xbar[n]:
								xbar[n][(i,j,t)] += r._value[1][n][(i,j,t)]
							else:
								xbar[n][(i,j,t)] = r._value[1][n][(i,j,t)]
				except:
					print(Tree_Graph)
					pdb.set_trace()
	### Ws needs to be calculated using xbar
	ws = {}
	for pd in results:
		for s in pd._value[0]:
			ws[s] = {}
			for n in pd._value[0][s]:
				for (i,j,t) in pd._value[0][s][n]:
					ws[s][(i,j,t)] = rho * (pd._value[0][s][n][(i,j,t)] - xbar[n][(i,j,t)])
			
	while alpha > tol and itr < max_iterations:
		
		### Create Pool
		pool = mp.Pool(ng)
		
		setlist = _grouper(scenarios,np)		
		results = [pool.apply_async(ss_PH_solve, args=(ss,rho,ws,xbar)) for ss in setlist]
		pool.close()
		pool.join()
		
		
		
		### Convert to old ws and xbar
		ws_old = copy.deepcopy(ws)
	
		### Xbar needs to be combine
		xbar = {}
		for n in builtins.Tree_Graph_GLOBAL.G.node:
			if n.startswith('s'):
				pass
			else:
				if n not in xbar:
					xbar[n] = {}
				for r in results:
					if n in r._value[1]:
						for (i,j,t) in r._value[1][n]:
							if (i,j,t) in xbar[n]:
								xbar[n][(i,j,t)] += r._value[1][n][(i,j,t)]
							else:
								xbar[n][(i,j,t)] = r._value[1][n][(i,j,t)]
		
		### Ws needs to be calculated using xbar
		ws = {}
		ENPV = 0
		a = 0
		for pd in results:
			ENPV += pd._value[2]
			for s in pd._value[0]:
				ws[s] = {}
				gk = 0
				for n in pd._value[0][s]:
					for (i,j,t) in pd._value[0][s][n]:
						ws[s][(i,j,t)] = ws_old[s][(i,j,t)] + rho * (pd._value[0][s][n][(i,j,t)] - xbar[n][(i,j,t)])
						gk += (pd._value[0][s][n][(i,j,t)] - xbar[n][(i,j,t)])**2
				gk = sqrt(gk)
				a += gk
		
		
		alpha = a
		itr += 1
	
	return ENPV

def ss_PH_solve(scenarios,rho,ws_o,xbar_o):
	ws = {}
	xbar = {}
	tENPV = 0 
	for s in scenarios:
		if s != None:
			
			### Calculate Node Names
			node_names = networkx.all_simple_paths(builtins.Tree_Graph_GLOBAL.G,'R',s)
			node_names = list(node_names)
			node_names = node_names[0]
			
			
			### Construct Model
			model,nodeset,nodevar,nodeallvar = model_constructor(s,node_names)
			
			model.obj.expr = model.obj.expr + sum(ws_o[s][(i,j,t)] * nodevar[n][str(n) + 'VarX'][i,j,t] for n in node_names for i in model.I for j in model.J for t in model.T if (i,j,t) in nodeset[n][str(n)+'SetX']) +\
								+ rho/2 * sum((nodevar[n][str(n) + 'VarX'][i,j,t] - xbar_o[n][(i,j,t)])**2 for n in node_names for i in model.I for j in model.J for t in model.T if (i,j,t) in nodeset[n][str(n)+'SetX'])
		
				
			model.preprocess()
			
			### Optimization things
			opt = SolverFactory("cplex")
			opt.options['threads'] = 3
			
			### Solve
			results = opt.solve(model)
			model.solutions.load_from(results)
			
			### Store Ws and xbar
			ws[s] = {}
			for n in node_names:
				
				if n.startswith('s'):
					pass
				else:
					if n not in ws[s]:
						ws[s][n] = {}
					if n not in xbar:
						xbar[n] = {}
						
					for (i,j,t) in nodeset[str(n)][str(n) +'SetX']:
						
						### Store Xbar
						if (i,j,t) in xbar[n]:
							
							### Calc Branch Probability
							Branch_PB = Calc_Branch_PB(n)
							
							### Update if it already exists
							xbar[n][(i,j,t)] += builtins.Tree_Graph_GLOBAL.G.node[s]['PB']/Branch_PB * nodevar[str(n)][str(n) + 'VarX'][i,j,t].value	
						
						else:
							### We have to add it if it hasn't been added before
							xbar[n][(i,j,t)] = 0
							
							### Calc Branch Probability
							Branch_PB = Calc_Branch_PB(n)
							
							
							xbar[n][(i,j,t)] += builtins.Tree_Graph_GLOBAL.G.node[s]['PB']/Branch_PB * nodevar[str(n)][str(n) + 'VarX'][i,j,t].value
							
						
						### Store ws
						ws[s][n][(i,j,t)] = nodevar[str(n)][str(n) + 'VarX'][i,j,t].value
			
			PB = builtins.Tree_Graph_GLOBAL.G.node[s]['PB']
		
			tENPV += PB * -model.obj()
			
	return ws, xbar, tENPV

def fs_PH_solve(scenarios,rho):
	
	ws = {}
	xbar = {}
	
	for s in scenarios:
		if s != None:
			
			### Calculate Node Names
			node_names = networkx.all_simple_paths(builtins.Tree_Graph_GLOBAL.G,'R',s)
			node_names = list(node_names)
			node_names = node_names[0]
			
			
			### Construct Model
			model,nodeset,nodevar,nodeallvar = model_constructor(s,node_names)
			
			### Optimization things
			opt = SolverFactory("cplex")
			opt.options['threads'] = 3
			
			### Solve
			results = opt.solve(model)
			model.solutions.load_from(results)
			
			
			### Store Ws and xbar
			ws[s] = {}
			for n in node_names:
				
				if n.startswith('s'):
					pass
				else:
					if n not in ws[s]:
						ws[s][n] = {}
					if n not in xbar:
						xbar[n] = {}
						
					for (i,j,t) in nodeset[str(n)][str(n) +'SetX']:
						
						### Store Xbar
						if (i,j,t) in xbar[n]:
							
							### Calc Branch Probability
							Branch_PB = Calc_Branch_PB(n)
							
							### Update if it already exists
							xbar[n][(i,j,t)] += builtins.Tree_Graph_GLOBAL.G.node[s]['PB']/Branch_PB * nodevar[str(n)][str(n) + 'VarX'][i,j,t].value	
						
						else:
							### We have to add it if it hasn't been added before
							xbar[n][(i,j,t)] = 0
							
							### Calc Branch Probability
							Branch_PB = Calc_Branch_PB(n)
							
							
							xbar[n][(i,j,t)] += builtins.Tree_Graph_GLOBAL.G.node[s]['PB']/Branch_PB * nodevar[str(n)][str(n) + 'VarX'][i,j,t].value
							
						
						### Store ws
						ws[s][n][(i,j,t)] = nodevar[str(n)][str(n) + 'VarX'][i,j,t].value
	
	
	return ws, xbar		
						
					
def Calc_Branch_PB(n):
	rpath = networkx.all_simple_paths(builtins.Tree_Graph_GLOBAL.G,'R',n)
	
	rpath = list(rpath)
	
	if len(rpath) > 1:
		print("Error in Path")
		exit()
	
	elif len(rpath) == 0:
		### This occurs when the node = R
		return 1
	
	else:
		path = rpath[0]
		nd = 0
		Branch_PB = 1
		while nd < len(path) - 1:
			Branch_PB = Branch_PB * builtins.Tree_Graph_GLOBAL.G.edge[path[nd]][path[nd+1]]['probability']
			nd += 1
		return Branch_PB
				

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

def model_constructor(scenario_name, node_names):
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
	
	return m, nodeset, nodevar, nodeallvar

def _evaluate_ENPV(model,node_names,nodeset,nodevar,nodeallvar,model_data,s):
	### Build Box
	tbox = [0 for _ in range(len(model_data._data['time_step'][None]))]
	jbox = [copy.deepcopy(tbox) for _ in range(len(model_data._data['trial'][None]))]
	xbox = [copy.deepcopy(jbox) for _ in range(len(model_data._data['product'][None]))]
	zbox = [copy.deepcopy(jbox) for _ in range(len(model_data._data['product'][None]))]		
	
	for n in node_names:
		for i in model_data._data['product'][None]:
			idx = model_data._data['product'][None].index(i)
			for j in model_data._data['trial'][None]:
				jdx = model_data._data['trial'][None].index(j)
				t = builtins.Tree_Graph_GLOBAL.G.node[n]['Active Periods']
				tdx = model_data._data['time_step'][None].index(t)
				### Determine where i,j,t is
				if n.startswith('s'):
					if (i,j,t) in nodeset[n][str(n) + 'SetXX']:
						xbox[idx][jdx][tdx] = round(nodevar[n][str(n)+'VarXX'][i,j,t].value,6)
						zbox[idx][jdx][tdx] = round(nodeallvar[n]['Z'][i,j,t].value,6)
				else:
					
					if (i,j,t) in nodeset[n][str(n) + 'SetX']:
						xbox[idx][jdx][tdx] = round(nodevar[n][str(n)+'VarX'][i,j,t].value,6)
						zbox[idx][jdx][tdx] = round(nodeallvar[n]['Z'][i,j,t].value,6)
					elif (i,j,t) in nodeset[n][str(n) + 'SetXX']:
						xbox[idx][jdx][tdx] = round(nodevar[n][str(n)+'VarXX'][i,j,t].value,6)
						zbox[idx][jdx][tdx] = round(nodeallvar[n]['Z'][i,j,t].value,6)
					else:
						xbox[idx][jdx][tdx] = return_fixed_values(i,j,t,n)
						zbox[idx][jdx][tdx] = round(nodeallvar[n]['Z'][i,j,t].value,6)
						
	ENPV = _Calculate_Value(xbox,zbox,model_data,s)
	
	return ENPV

def _Calculate_Value(X, Z,model_data,scenario_name):
	### Use M2S Object to calculate parameters
	
	import Core.Solvers.MTSSP.M2S_item as M2S_item
	duration = model_data._data['trial_duration']
	gammaL = {}
	gammaD = {}
	revenue_max = {}
	prod = model_data._data['product'][None]
	sg = model_data._data['trial'][None]
	ts = model_data._data['time_step'][None]
	resource_type = model_data._data['resource_type'][None]
	
	trial_cost = model_data._data['trial_cost']
	
	success = {}
	for i in builtins.MD_GLOBAL._data['product'][None]:
		if builtins.Tree_Graph_GLOBAL.G.node[scenario_name]['Outcome'][builtins.MD_GLOBAL._data['product'][None].index(i)] == len(sg):
			success[i] = 1
		else:
			success[i] = 0
	
	for items in model_data._data['gammaL']:
		gammaL[items[0]] = model_data._data['gammaL'][items]
		
	for items in model_data._data['gammaD']:
		gammaD[items[0]] = model_data._data['gammaD'][items]

	for items in model_data._data['maximum_revenue']:
		revenue_max[items[0]] = model_data._data['maximum_revenue'][items]
		
	last_trial = sg[-1]
	
	## Calculate running rev
	rev_run = M2S_item.calc_rr(revenue_max,gammaL,duration, prod, sg, ts)

	##Calculate open rev
	rev_open = M2S_item.calc_openrev(revenue_max,gammaL,duration, prod, sg, ts, len(ts))

	##Calculate Discounting Factor
	discounting_factor = M2S_item.calc_discounting_factor(revenue_max,gammaL,trial_cost, prod, sg, len(ts))
	
	
	cost = sum((1 - 0.025 * (t - 1)) * trial_cost[(i,j)]*X[prod.index(i)][sg.index(j)][ts.index(t)] for i in prod for j in sg for t in ts)
	
	revenue =  sum( success[i] * revenue_max[i]*X[prod.index(i)][sg.index(last_trial)][ts.index(t)] for i in prod for t in ts) -  sum(success[i]*gammaD[i]*Z[prod.index(i)][sg.index(j)][ts.index(t)] for i in prod for j in sg if j>1 for t in ts) - sum(gammaL[i]*success[i]*(t + duration[(i,len(sg))])*X[prod.index(i)][sg.index(last_trial)][ts.index(t)] for i in prod for t in ts)
	
	FRV1 = sum(success[i] * rev_open[(i,j)] * discounting_factor[(i,j)] * Z[prod.index(i)][sg.index(j)][ts.index(len(ts))] for i in prod for j in sg)

	FRV2 = sum(success[i] * rev_run[(i,j,t)] * discounting_factor[(i,j+1)] * X[prod.index(i)][sg.index(j)][ts.index(t)] for i in prod for j in sg if j < last_trial for t in ts if t > len(ts) - duration[(i,j)])

	Free_revenue = FRV1 + FRV2
	
	
	value = revenue + Free_revenue - cost
	
	return value
	
def _grouper(ibl, n, fillvalue=None):
	args = [iter(ibl)]*n
	return itertools.zip_longest(fillvalue = fillvalue, *args)
