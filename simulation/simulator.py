import sys
import math
from os import path
sys.path.append( path.dirname(path.dirname( path.abspath(__file__)) ) ) 
import settings

def simulate(system, mode='section',devMultiplier=1,start='first',steps='all',flow='b_b',VLE='PR'):
	"""
	modes can be section or devision.
	devMultiplier is the size of step for running simulation in case of devision
	steps can be "all" or an int
	start can be "last" or "first"
	"""
	assert mode=="section" or mode=="devision","mode must be either 'section' or 'devision'"
	assert isinstance(system,settings.System),"system must be an instance of the System class"
	assert type(devMultiplier)==int,"devMultiplier must be an integar"
	assert start=="last" or start=="first","start must be either 'last' or 'first'"
	assert steps=="all" or type(steps)==int,"steps must be either 'all' or an int"
	simData=system.simData
	# {
	# 	'vars':[
			# {
			# 	'geo_index':None,'T':None,'P':None,'v_g':None,'v_l':None,'vapQ':None,'r_g':None,'r_l':None,'m_g':None,'m_l':None,'pattern':None,'hold_up':None,'Cp_l':None,'Cp_g':None
			# }
	# 	],'geometry':[
	# 		{'index':None,'x':None,'y':None,'d':None}
	# 	]
	# }
	if start=='first':
		position=0
		if system.pvtType=='bo':
			simData['vars']=[{
				'geo_index':0,
				'T':system.T_in,
				'P':system.p_in,
				'v_g':None,
				'v_l':None,
				'vapQ':system.pvt['vapQ'],
				'r_g':system.pvt['r_g'],
				'r_l':system.pvt['r_l'],
				'm_g':system.pvt['m_g'],
				'm_l':system.pvt['m_l'],
				'pattern':None,
				'hold_up':None,
				'Cp_l':None,
				'Cp_g':None
				}]
		elif system.pvtType=='comp':
			#calculate initial pvt params
			pass
		else:
			assert False, 'pvt type error'
	elif start=='last':
		position=len(simData['vars']['P'])
		if position==0:
			if system.pvtType=='bo':
				simData['vars']=[{
				'geo_index':0,
				'T':system.T_in,
				'P':system.p_in,
				'v_g':None,
				'v_l':None,
				'vapQ':system.pvt['vapQ'],
				'r_g':system.pvt['r_g'],
				'r_l':system.pvt['r_l'],
				'm_g':system.pvt['m_g'],
				'm_l':system.pvt['m_l'],
				'pattern':None,
				'hold_up':None,
				'Cp_l':None,
				'Cp_g':None
				}]
			elif system.pvtType=='comp':
				#calculate initial pvt params
				pass
	else:
		assert False
	stepIndices=[]
	if steps=='all':
		end=len(simData['geometry'])-1
		stepIndices.append(position)
		if mode=='section':
			end=position
			while end<len(simData['geometry']):
				if not simData['geometry'][end]['index']-simData['geometry'][end-1]['index']==0:
					if not end==position:
						stepIndices.append(end)
				end+=1
			if simData['geometry'][end-1]['index']-simData['geometry'][end-2]['index']==0:
				stepIndices.append(end-1)
		elif mode=='devision':
			end=position
			while end<len(simData['geometry']):
				if not end==position:
					stepIndices.append(end)
				end+=devMultiplier
			if end>=len(simData['geometry']) and not end-devMultiplier==len(simData['geometry'])-1:
				stepIndices.append(len(simData['geometry'])-1)
		else:
			assert False

	else:
		if mode=='section':
			end=position
			stepIndices.append(position)
			s=steps
			while end<len(simData['geometry']) and s>0:
				if not simData['geometry'][end]['index']-simData['geometry'][end-1]['index']==0:
					if not end==position:
						s=s-1
						stepIndices.append(end)
				end+=1
			if end==len(simData['geometry']):
				stepIndices.append(end-1)
		elif mode=='devision':
			end=position
			stepIndices.append(position)
			s=steps
			while end<len(simData['geometry']) and s>0:
				if not end==position:
					s=s-1
					stepIndices.append(end)
				end+=devMultiplier
			if end>=len(simData['geometry']) and not end-devMultiplier==len(simData['geometry'])-1:
				stepIndices.append(len(simData['geometry'])-1)
		else:
			assert False

	print stepIndices