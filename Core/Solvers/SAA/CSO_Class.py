import os
import sys
import pdb

class CSO:
	def __init__(self,name, lists_of_CLists, Scenarios):
		self.Name = name
		
		self.CLists = lists_of_CLists
		
		#### Create Set List
		self.sets = self._create_setlists(Scenarios)
		
	def __repr__(self):
		return self.name
		
	def _create_setlists(self,Scenarios):
		sets = {}
		for j in self.CLists:
			sets[tuple(j)] = [[i,] for i in Scenarios]
		
		return sets
