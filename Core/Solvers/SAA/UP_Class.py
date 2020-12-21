import os
import sys

class Uncertain_Parameter:
	
	def __init__(self, Name, Type, Realization_Type, Realizations, Realization_Order = [], GRealization_Order_Sets = []):
		### This is the name of the uncertain parameter
		self.Name = Name
		
		### Uncertain parameters are either endogenous or exogenous
		self.Type = Type
		
		### Realizations of the uncertain parameter are either gradual or instant
		self.Realization_Type = Realization_Type
		
		### Values for the realizations
		self.Realizations = Realizations
		
		### Uncertain parameters that must be realized before this parameter can be realized
		self.ROrder = Realization_Order 
		
		### Gradual realization order and sets (i.e. [0,1,[2,3]])
		self.GROrderSets = GRealization_Order_Sets

	def __repr__(self):
		return self.Name
