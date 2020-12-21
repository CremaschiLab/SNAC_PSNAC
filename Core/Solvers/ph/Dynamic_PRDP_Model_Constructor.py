### Progressive hedging requires a model with specified first and second stage decisions.
### The premise of this constructor is to construct a model where first and second stage decisions are specified

import sys
import os
from os.path import abspath, dirname
sys.path.insert(0, dirname(dirname(dirname(dirname(abspath(__file__))))))
from pyomo.environ import *
from pyomo.dae import *
import pdb

def construct_model(sets,params,success, Fixed_Variables, Constrained_Variables):

	m = ConcreteModel()
	
	##Sets##
	m.I = Set(initialize=sets['Product'])
	m.J = Set(initialize=sets['Trial'])
	m.T = Set(initialize=sets['Time_Step'])
	m.R = Set(initialize=sets['resource_type'])
	
	def FS_IDX_init(model):
		return ((i,j,t) for i in m.I for j in m.J for t in m.T if (i,j,t) in Constrained_Variables)
	m.FS_IDX = set(within = m.I * m.J * m.T, initialize = FS_IDX_init)
	
	def SS_IDX_init(model):
		return ((i,j,t) for i in m.I for j in m.J for t in m.T if (i,j,t) not in Constrained_Variables and not in Fixed_Variables)
	m.SS_IDX = set(within = m.I * m.J * m.T, initialize = SS_IDX_init)
	
	
	### Fixed Parameters###
	m.Resource_Max = Param(initialize=params['resource_max'])
	m.GammaL = Param(initialize = params['gammaL'])
	m.GammaD = Param(initialize = params['gammaD'])
	m.Duration = Param(initialize = params['duration']) 
	m.Trial_Cost = Param(initialize = params['trial_cost']) 
	m.Resources_Required = Param(initialize = params['resource_required'])
	m.Revenue_Max = Param(initialize = params['revenue_max'])
	m.LastTimeStep = Param(initialize = params['last_time_step'])
	m.Last_Trial = Param(initialize = params['last_trial'])
	m.RR = Param(initialize = params['Running_Revenue'])
	m.OR = Param(initialize = params['Open_Revenue'])
	m.DF = Param(initialize = params['Discounting_Factor'])
	
	### Scenario Specific Parameters###
	m.Success = Param(initialize = success, muteable = True)
	
	### Variables ###
	m.FS_X = Var(m.FS_IDX, bounds=(0,1))
	m.SS_X = Var(m.SS_IDX, bounds=(0,1))
	m.Y = Var(m.I,m.J,m.T, bounds=(0,1))
	m.Z = Var(m.I,m.J,m.T, bounds= (0,1))
	
	##############################################################################
	### Constraints
	##############################################################################
	
	
	### Constraint-- The trial is complete if the trial was started a duration ago
	def Trial_Finish_rule(m,i,j,t):
		past_duration = t - m.Duration[i,j]
		if t==1:
			return model.Decision_Y[i,j,t] == 0
		elif  t > 1 and  past_duration < 1 :
			return model.Decision_Y[i,j,t] == model.Decision_Y[i,j,t-1]
		else:
			if (i,j,t) in Constrained_Variables:
				return m.Y[i,j,t] == m.Y[i,j,t-1] + m.FS_X[i,j,past_duration] 
			else:
				return m.Y[i,j,t] == m.Y[i,j,t-1] + m.SS_X[i,j,past_duration] 
	m.Trial_Finish = Constraint(m.I, m.J, m.T, rule=Trial_Finish_rule)

	def Calculate_Z_rule(m,i,j,t):
		if (i,j,t) in Fixed_Variables:
			decX = Fixed_Variables[i,j,t]
			if j > 1:
				pd = t - m.duration[i,j-1,t]
				if pd > 0 and (i,j-1,pd) in Fixed_Variables:
					decX_pd = Fixed_Variables[(i,j,pd)]
				elif pd > 0 and (i,j-1,pd) in Constrained_Variables:
					decX_pd = m.FS_X[i,j,pd]
				elif pd > 0:
					decX_pd = m.SS_X[i,j,pd]
				
		elif (i,j,t) in Constrained_Variables
			decX = m.FS_X[(i,j,t)]
			if j > 1:
				pd = t - m.duration[i,j-1,t]
				if pd > 0 and (i,j-1,pd) in Fixed_Variables:
					decX_pd = Fixed_Variables[(i,j,pd)]
				elif pd > 0 and (i,j-1,pd) in Constrained_Variables:
					decX_pd = m.FS_X[i,j,pd]
				elif pd > 0:
					decX_pd = m.SS_X[i,j,pd]
		else:
			decX = m.SS_X[(i,j,t)]
			if j > 1:
				pd = t - m.duration[i,j-1,t]
				if pd > 0 and (i,j-1,pd) in Fixed_Variables:
					decX_pd = Fixed_Variables[(i,j,pd)]
				elif pd > 0 and (i,j-1,pd) in Constrained_Variables:
					decX_pd = m.FS_X[i,j,pd]
				elif pd > 0:
					decX_pd = m.SS_X[i,j,pd]
			
		if j == 1:
			if t > 1:
				return m.Z[i,j,t] == m.Z[i,1,t-1] - decX
			else:
				return m.Z[i,j,1] == 1 - decX
		else:
			if t- m.Duration[i,j-1] >0 and t>1:
				return m.Z[i,j,t] == m.Z[i,j,t-1] + decX_pd - decX
			elif t-m.Duration[i,j-1] <= 0 and t > 1:
				return m.Z[i,j,t] == m.Z[i,j,t-1] - decX
			elif t==1:
				return m.Z[i,j,t] == -decX
	m.Calculate_Z = Constraint(m.I,m.J,m.T, rule= Calculate_Z_rule) 

	### Constraint--You can only start each trial once
	def Constraint_3_rule(model,i,j):
		exp = 0
		for t in m.T:
			if (i,j,t) in Fixed_Variables:
				exp += Fixed_Variables[i,j,t]
			elif (i,j,t) in Constrained_Variables:
				exp += m.FS_X[i,j,t]
			else:
				exp += m.SS_X[i,j,t]
 		return exp <= 1
	model.Constraint_3 = Constraint(m.I,m.J,rule=Constraint_3_rule)

	### Constraint-- Must Start Previous Trial Before it can be finished
	def Constraint_4_rule(model,i,j,t):
		if j > 1:
			previous_trial = j-1
			LHS = 0
			for tprime in m.T:
				if tprime =< t:
					if (i,j,tprime) in Fixed_Variables:
						LHS += Fixed_Variables[(i,j,tprime)] 
					elif (i,j,tprime) in Constrained_Variables:
						LHS += m.FS_X[i,j,tprime]
					else:
						LHS += m.SS_X[i,j,tprime]
			return LHS <= model.Decision_Y[i,previous_trial,t]
		else:
			return Constraint.Skip
	model.Constraint_4 = Constraint(m.I, m.J, m.T, rule=Constraint_4_rule)

	### Constraint-- Ensures resources are managed correctly
	def Resource_Constraint_rule(model,r,t):
		LHS = 0
		for i in m.I:
			for j in m.J:
				for tprime in m.T:
					if tprime > (t - m.Duration[i,j]) and tprime <= t:
						if (i,j,tprime) in Fixed_Variables:
							LHS += Fixed_Variables[i,j,tprime] * m.RR[i,j,r]
						elif (i,j,tprime) in Constrained_Variables:
							LHS += m.FS_X[i,j,tprime]
						else:
							LHS += m.SS_X[i,j,tprime]
		return LHS <= model.Resource_Max[r]
	model.Resource_Constraint = Constraint(m.R, m.T, rule= Resource_Constraint_rule)
	
	def ComputeFirstStageCosts_rule():
		cst = 0
		### Costs of Starting First Stage Trials
		cst += -sum((1 - (0.025 * (t - 1))) * m.Trial_Cost[i,j] * m.FS_X[i,j,t] for (i,j,t) in m.FS_IDX)
		
		cst += -sum((1 - (0.025 * (t - 1))) * m.Trial_Cost[i,j] * Fixed_Variables[(i,j,t)] for (i,j,t) in Fixed_Variables)
		
		### Revenue From Starting Final Trials
		cst += sum(m.Success[i]*(m.Revenue_Max[i] * m.FS_X[i,j,t]) for (i,j,t) in m.FS_IDX if j == m.last_trial)
		
		cst += sum(m.Success[i]*(m.Revenue_Max[i] * Fixed_Variables[(i,j,t)]) for (i,j,t) in Fixed_Variables if j == m.last_trial)
		
		### Losses From Not Starting Trials
		cst += sum(-model.Success[i]*model.GammaL[i]*(t + m.Duration[i, j]) * m.FS_X[i,j,t] for (i,j,t) in m.FS_IDX if j == m.last_trial)
		
		cst += sum(-model.Success[i]*model.GammaL[i]*(t + m.Duration[i, j]) * Fixed_Variables[(i,j,t)] for (i,j,t) in Fixed_Variables if j == m.last_trial)
		
		### Free Revenue From Products Started But Not Yet Finished
		cst	+= sum(m.Success[i] * m.RR[i,j,t] * m.DF[i,j+1]*m.FS_X[i,j,t] for (i,j,t) in m.FS_IDX if j < m.Last_Trial for t in m.T if t > m.LastTimeStep - m.Duration[i,j])
		
		cst	+= sum(m.Success[i] * m.RR[i,j,t] * m.DF[i,j+1]*Fixed_Variables[(i,j,t)] for (i,j,t) in Fixed_Variables if j < m.Last_Trial for t in m.T if t > m.LastTimeStep - m.Duration[i,j])
		
		return -cst
	model.FirstStageCost = Expression(rule=ComputeFirstStageCost_rule)
	
	def ComputeSecondStageCost_rule():
		cst = 0
		
		### Costs of Starting First Stage Trials
		cst += -sum((1 - (0.025 * (t - 1))) * model.Trial_Cost[i,j] * m.SS_X[i,j,t] for (i,j,t) in m.SS_IDX)
	
		### Revenue From Starting Final Trials
		cst += sum(m.Success[i]*(m.Revenue_Max[i] * m.SS_X[i,j,t]) for (i,j,t) in m.SS_IDX if j == m.last_trial)
		
		### Losses For Decisions Not Made
		cst += sum(-m.Success[i]*m.GammaD[i] * sum(m.Z[i,j,t] for j in m.J if j > 1) for i in m.I for t in m.T)
	
		### Losses for the time to development
		cst += sum(-m.Success[i]*m.GammaL[i]*(t + m.Duration[i, j]) * m.SS_X[i,j,t] for (i,j,t) in m.SS_IDX if j == m.Last_Trial)
	
		### Revenue Yet to Be Realized
		cst += sum(m.Success[i] * m.OR[i,j] * m.DF[i,j]* m.Z[i,j,m.LastTimeStep] for i in m.I for j in m.J)
		
		### Revenue for Started Products Not Yet Realized
		cst += sum(m.Success[i] * m.RR[i,j,t] * m.DF[i,j+1]*m.SS_X[i,j,t] for (i,j,t) in m.SS_IDX if j < m.Last_Trial for t in m.T if t > m.LastTimeStep - m.Duration[i,j])
		
		return -cst
	model.SecondStageCost = Expression(rule=ComputeSecondStageCost_rule)
	
	return model
	
def multistage_constructor(Tree_Object, FV, CD):
	
	m = ConcreteModel()
	
	##Sets##
	m.I = Set(initialize=sets['Product'])
	m.J = Set(initialize=sets['Trial'])
	m.T = Set(initialize=sets['Time_Step'])
	m.R = Set(initialize=sets['resource_type'])
	
	def IDX_init(model,indicies):
		return ((i,j,t) for i in m.I for j in m.J for t in m.T if (i,j,t) in indicies)
		
	for node_name in Tree_Object.nodes:
		m.add_component(node_name, Set(within = m.I * m.J * m.T, initialize = IDX_init(node_name))
		
	### Fixed Parameters###
	m.Resource_Max = Param(initialize=params['resource_max'])
	m.GammaL = Param(initialize = params['gammaL'])
	m.GammaD = Param(initialize = params['gammaD'])
	m.Duration = Param(initialize = params['duration']) 
	m.Trial_Cost = Param(initialize = params['trial_cost']) 
	m.Resources_Required = Param(initialize = params['resource_required'])
	m.Revenue_Max = Param(initialize = params['revenue_max'])
	m.LastTimeStep = Param(initialize = params['last_time_step'])
	m.Last_Trial = Param(initialize = params['last_trial'])
	m.RR = Param(initialize = params['Running_Revenue'])
	m.OR = Param(initialize = params['Open_Revenue'])
	m.DF = Param(initialize = params['Discounting_Factor'])
	
	### Scenario Specific Parameters###
	m.Success = Param(initialize = success, muteable = True)

	### Variables ###
	m.X = Var(m.FS_IDX, bounds=(0,1))
	m.XX = Var(m.SS_IDX, bounds=(0,1))
	mvars = [m.X, m.XX]
	m.Y = Var(m.I,m.J,m.T, bounds=(0,1))
	m.Z = Var(m.I,m.J,m.T, bounds= (0,1))
	
	
	def ComputeStageCost_rule():
		
		cst = 0
		
		for i in range(2):
			decvar = mvars[i]
			
			### Costs of Starting First Stage Trials
			cst += -sum((1 - (0.025 * (t - 1))) * model.Trial_Cost[i,j] * decvar[i,j,t] for (i,j,t) in m.SS_IDX ### this should be the stageX indicies)
		
			### Revenue From Starting Final Trials
			cst += sum(m.Success[i]*(m.Revenue_Max[i] * decvar[i,j,t]) for (i,j,t) in m.SS_IDX if j == m.last_trial)
			
			### Losses for the time to development
			cst += sum(-m.Success[i]*m.GammaL[i]*(t + m.Duration[i, j]) * decvar[i,j,t] for (i,j,t) in m.SS_IDX if j == m.Last_Trial)
		
			### Revenue for Started Products Not Yet Realized
			cst += sum(m.Success[i] * m.RR[i,j,t] * m.DF[i,j+1]*decvar[i,j,t] for (i,j,t) in m.SS_IDX if j < m.Last_Trial for t in m.T if t > m.LastTimeStep - m.Duration[i,j])
		
		
		### Losses For Decisions Not Made
		cst += sum(-m.Success[i]*m.GammaD[i] * sum(m.Z[i,j,t] for j in m.J if j > 1) for i in m.I for t in m.T)
			
		### Revenue Yet to Be Realized
		cst += sum(m.Success[i] * m.OR[i,j] * m.DF[i,j]* m.Z[i,j,m.LastTimeStep] for i in m.I for j in m.J)
		
		return -cst
	

	
