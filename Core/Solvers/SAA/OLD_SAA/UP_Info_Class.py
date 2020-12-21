import os
import sys

class Uncertain_Parameter:
	
	def __init__(self, Name, Type, Realization_Type, Realizations):
		Uncertain_Parameter.Name = Name
		Uncertain_Parameter.Type = Type
		Uncertain_Parameter.Realization_Type = Realization_Type
		Uncertain_Parameter.Realizations = Realizations

	def __repr__(self):
		return self.Name
