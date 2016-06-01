import sys
import math
from os import path
sys.path.append( path.dirname(path.dirname( path.abspath(__file__)) ) ) 
import settings

def simulate(system, mode='section',devMultiplier=1,start='first',steps='all'):
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
	#{'t':[],'x':[],'d':[],'y':[],'i':[],'vars':{'T':[],'P':[],'v_g':[],v_l,'vapQ':[],'r_g':[],'r_l':[],'m_g':[],'m_l':[]}}
	if start=='first':
		position=0
		if system.pvtType=='bo':
			Ql=(1-system.pvt['vapQ'])*system.massFlow/system.pvt['r_l']
			Qg=system.pvt['vapQ']*system.massFlow/system.pvt['r_g']
			H=Ql/(Ql+Qg) ####No Slip assumed
			print H
			ID=system.sections[simData['i'][position]].d-2*system.sections[simData['i'][position]].t
			area=math.pi*math.pow(ID,2)/4
			simData['T']=[system.T_in,]
			simData['P']=[system.p_in,]
			simData['v_g']=[Qg/((1-H)*area),]
			simData['v_l']=[Ql/(H*area),]
			simData['vapQ']=[system.pvt['vapQ'],]
			simData['r_l']=[system.pvt['r_l'],]
			simData['r_g']=[system.pvt['r_g'],]
			simData['m_g']=[system.pvt['m_g'],]
			simData['m_l']=[system.pvt['m_l'],]
		elif system.pvtType=='comp':
			#calculate initial pvt params
			pass
	elif start=='last':
		position=len(simData['vars']['P'])
		if position==0:
			if system.pvtType=='bo':
				Ql=(1-system.pvt['vapQ'])*system.massFlow/system.pvt['r_l']
				Qg=system.pvt['vapQ']*system.massFlow/system.pvt['r_g']
				H=Ql/(Ql+Qg) ####No Slip assumed
				ID=system.sections[simData['i'][position]].d-2*system.sections[simData['i'][position]].t
				area=math.pi*math.pow(ID,2)/4
				simData['T']=[system.T_in,]
				simData['P']=[system.p_in,]
				simData['v_g']=[Qg/((1-H)*area),]
				simData['v_l']=[Ql/(H*area),]
				simData['vapQ']=[system.pvt['vapQ'],]
				simData['r_l']=[system.pvt['r_l'],]
				simData['r_g']=[system.pvt['r_g'],]
				simData['m_g']=[system.pvt['m_g'],]
				simData['m_l']=[system.pvt['m_l'],]
			elif system.pvtType=='comp':
				#calculate initial pvt params
				pass
	else:
		assert False
	stepIndices=[]
	if steps=='all':
		end=len(simData['d'])-1
		stepIndices.append(position)
		if mode=='section':
			end=position
			while end<len(simData['d']):
				if not simData['i'][end]-simData['i'][end-1]==0:
					if not end==position:
						stepIndices.append(end)
				end+=1
			if simData['i'][end-1]-simData['i'][end-2]==0:
				stepIndices.append(end-1)
		elif mode=='devision':
			end=position
			while end<len(simData['d']):
				if not end==position:
					stepIndices.append(end)
				end+=devMultiplier
			if end>=len(simData['d']) and not end-devMultiplier==len(simData['d'])-1:
				stepIndices.append(len(simData['d'])-1)
		else:
			assert False

	else:
		if mode=='section':
			end=position
			stepIndices.append(position)
			s=steps
			while end<len(simData['d']) and s>0:
				if not simData['i'][end]-simData['i'][end-1]==0:
					if not end==position:
						s=s-1
						stepIndices.append(end)
				end+=1
			if end==len(simData['d']):
				stepIndices.append(end-1)
		elif mode=='devision':
			end=position
			stepIndices.append(position)
			s=steps
			while end<len(simData['d']) and s>0:
				if not end==position:
					s=s-1
					stepIndices.append(end)
				end+=devMultiplier
			if end>=len(simData['d']) and not end-devMultiplier==len(simData['d'])-1:
				stepIndices.append(len(simData['d'])-1)
		else:
			assert False

		print stepIndices