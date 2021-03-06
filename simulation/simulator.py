import math
import settings
from simulation.fluid_flow import b_and_b as b_b
from simulation.fluid_flow import p_and_a as p_a
from simulation.fluid_flow import one_phase as o_p
# from simulation.thermodynamics import glaso
from simulation.thermodynamics import PVT
from simulation.heat_transfer import ambient_heat_transfer as ambient_heat_transfer
def calcStepIndices(geometry,devs):
	stepIndices=[]
	devisions=len(geometry)
	devMultiplier=devisions/devs
	mode='section' if 10>devs else 'devision'
	position=0
	end=len(geometry)-1
	stepIndices.append(position)
	if mode=='section':
		end=position
		while end<len(geometry):
			if not geometry[end]['index']-geometry[end-1]['index']==0:
				if not end==position:
					stepIndices.append(end)
			end+=1
		if geometry[end-1]['index']-geometry[end-2]['index']==0:
			stepIndices.append(end-1)
	elif mode=='devision':
		end=position
		while end<len(geometry):
			if not end==position:
				stepIndices.append(end)
			end+=devMultiplier
		if end>=len(geometry) and not end-devMultiplier==len(geometry)-1:
			stepIndices.append(len(geometry)-1)
	else:
		assert False
	return stepIndices

def forward_steady_itterate(system, tol_P=10,flow='b_b',VLE='ideal',devs=8, tol_T=2):
	#tol_P in Pascals
	assert type(tol_P)==int or type(tol_P)==float,"pressure convergence must be a number"
	assert isinstance(system,settings.System),"system must be an instance of the System class"
	assert flow in ['b_b','p_a']
	assert VLE in ['incompressible','ideal','pr']
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
	position=0
	if system.pvtType=='bo' and VLE=='incompressible':
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
			'Cp_l':system.pvt['Cp_l'],
			'Cp_g':system.pvt['Cp_g'],
			'enthalpy':None
			}]
	elif system.pvtType=='bo' and VLE=='ideal':
		simData['vars']=[{
			'geo_index':0,
			'T':system.T_in,
			'P':system.p_in,
			'v_g':None,
			'v_l':None,
			'vapQ':system.pvt['vapQ'],
			'r_g':system.pvt['r_g'],
			'r_l':system.pvt['r_l'],
			'm_g':1e-8,
			'm_l':system.pvt['m_l'],
			'pattern':None,
			'hold_up':None,
			'Cp_l':system.pvt['Cp_l'],
			'Cp_g':None,
			'enthalpy':None
			}]
	elif system.pvtType=='comp' and VLE=='pr':
		simData['vars'] = [PVT.get_props(system.pvt,system.T_in,system.p_in)]
		simData['vars'][0]['geo_index']= 0
		simData['vars'][0]['T']=system.T_in
		simData['vars'][0]['P']=system.p_in
		simData['vars'][0]['v_g']=None,
		simData['vars'][0]['v_l']=None,
		simData['vars'][0]['pattern']=None,
		simData['vars'][0]['hold_up']=None
	else:
		assert False, 'pvt type error'

	stepIndices=calcStepIndices(simData['geometry'],devs)
	geometry=simData['geometry']
	sim_vars=simData['vars']
	counter=0
	print sim_vars[0]
	for i in stepIndices:
		if counter==0:
			pass
		else:
			prev=sim_vars[-1]
			geo0=prev['geo_index']
			geo1=i
			z=geometry[geo1]['y']-geometry[geo0]['y']
			pipe_index=geometry[geo0]['index']
			pipe=system.sections[pipe_index]
			roughness=pipe.e
			thickness=pipe.t
			ID=pipe.d-2*thickness
			L=geometry[geo1]['d']-geometry[geo0]['d']
			U=pipe.U
			h0=prev['enthalpy']
			T0=prev['T']
			r_g0=prev['r_g']
			r_l0=prev['r_l']
			m_g0=prev['m_g']
			m_l0=prev['m_l']
			Cp_g0=prev['Cp_g']
			Cp_l0=prev['Cp_l']
			vapQ0=prev['vapQ']
			massFlow=system.massFlow
			P0=prev['P']
			print '****************************',vapQ0
			#pressure drop calc
			if vapQ0 < 1e-3:
				P1=o_p.calcPressureDrop(P0,massFlow,r_l0,m_l0,z,roughness,ID,L,T0)
				pattern = 'liquid'
				hold_up=1.
			elif vapQ0 > 1 - 1e-3:
				P1=o_p.calcPressureDrop(P0,massFlow,r_g0,m_g0,z,roughness,ID,L,T0)
				pattern='gas'
				hold_up=0.
			elif flow=='b_b':
				P1,hold_up,pattern=b_b.calcPressureDrop(P0,massFlow,vapQ0,r_g0,r_l0,m_g0,m_l0,z,roughness,ID,L,T0)
			elif flow=='p_a':
				try:
					P1,hold_up,pattern=p_a.calcPressureDrop(P0,massFlow,vapQ0,r_g0,r_l0,m_g0,m_l0,z,roughness,ID,L,T0)
				except:
					print  '<<<<<<<<<<<<<<<<<<<<<'
					print r_g0,r_l0
					P1, hold_up, pattern = b_b.calcPressureDrop(P0, massFlow, vapQ0, r_g0, r_l0, m_g0, m_l0, z,
																roughness, ID, L, T0)
			else:
				assert False,'invalid pressure drop method'
			# print 'P0: %s, P1: %s, T0: %s, L: %s, liq visc: %s, gas visc: %s, vapQ: %s, liq density: %s, gas density: %s' %(P0,P1,T0,L,m_l0,m_g0,vapQ0,r_l0,r_g0)
			assert P1>1, "pressure drop is too high"
			#temp drop calc
			T_out=system.T_amb
			T1=ambient_heat_transfer.calculate_T(ID+2*thickness,L,massFlow,Cp_l0,Cp_g0,vapQ0,T0,T_out,U)
			#pvt calc
			if VLE=='incompressible' and system.pvtType=='bo':
				r_g1=r_g0
				r_l1=r_l0
				m_g1=m_g0
				h1=None
				m_l1=m_l0 #water/air experiment
				#m_l1=glaso.updateLiquidVisc(m_l0,r_l1,T1,T0)
				vapQ1=vapQ0
				Cp_g1=Cp_g0
				Cp_l1=Cp_l0
			elif VLE=='ideal' and system.pvtType=='bo':
				r_g1=r_g0*P1*(T1+273.14)/(P0*(T0+273.14))
				r_l1=r_l0
				m_g1=m_g0
				h1=None
				m_l1=m_l0 #water/air experiment
				#m_l1=glaso.updateLiquidVisc(m_l0,r_l1,T1,T0)
				vapQ1=vapQ0
				Cp_g1=Cp_g0
				Cp_l1=Cp_l0
			elif VLE=='pr' and system.pvtType=='comp':
				props=PVT.get_props(system.pvt,T1,P1)
				r_g1=props['r_g']
				r_l1=props['r_l']
				m_g1=props['m_g']
				m_l1=props['m_l']
				vapQ1=props['vapQ']
				Cp_g1=props['Cp_g']
				Cp_l1=props['Cp_l']
				h1=props['enthalpy']

				# T_new=ambient_heat_transfer.minimize_q(ID+2*thickness,L,massFlow,.5*(Cp_l0+Cp_l1),0.5*(Cp_g1+Cp_g0),0.5*(vapQ0+vapQ1),T0,T1,T_out,U,h1,h0)
				T_new = ambient_heat_transfer.calculate_T(ID + 2 * thickness, L, massFlow, .5 * (Cp_l0 + Cp_l1),
														  0.5 * (Cp_g1 + Cp_g0), 0.5 * (vapQ0 + vapQ1), T0, T_out, U)

				T_error=T_new-T1
				T1=T_new
				while T_error>tol_T:
					props = PVT.get_props(system.pvt, T1, P1)
					r_g1 = props['r_g']
					r_l1 = props['r_l']
					m_g1 = props['m_g']
					m_l1 = props['m_l']
					vapQ1 = props['vapQ']
					Cp_g1 = props['Cp_g']
					Cp_l1 = props['Cp_l']
					h1 = props['enthalpy']

					# T_new=ambient_heat_transfer.minimize_q(ID+2*thickness,L,massFlow,.5*(Cp_l0+Cp_l1),0.5*(Cp_g1+Cp_g0),0.5*(vapQ0+vapQ1),T0,T1,T_out,U,h1,h0)
					T_new=ambient_heat_transfer.calculate_T(ID+2*thickness,L,massFlow,.5*(Cp_l0+Cp_l1),0.5*(Cp_g1+Cp_g0),0.5*(vapQ0+vapQ1),T0,T_out,U)
					T_error = T_new - T1

					T1=T_new
			prev['v_g'] = vapQ0 * massFlow * 4 / (ID * ID * math.pi * r_g0 * (1 - hold_up)) if hold_up < 1 else 0
			prev['v_l'] = (1 - vapQ0) * massFlow * 4 / (ID * ID * math.pi * r_l0 * hold_up) if hold_up > 0 else 0
			if flow in ['p_a','b_b']:
				try:
					prev2 = sim_vars[-2]
					print prev2
				except:
					print 'Fail'
					pass
				else:
					P1-=(r_g1*prev['v_g']**2-r_g0*prev2['v_g']**2)*.5
			prev['pattern']=pattern
			prev['hold_up']=hold_up
			sim_vars[-1]=prev
			sim_vars.append({
				'geo_index':geo1,
				'T':T1,
				'P':P1,
				'v_g':None,
				'v_l':None,
				'vapQ':vapQ1,
				'enthalpy':h1,
				'r_g':r_g1,
				'r_l':r_l1,
				'm_g':m_g1,
				'm_l':m_l1,
				'pattern':None,
				'hold_up':None,
				'Cp_l':Cp_l1,
				'Cp_g':Cp_g1
				})
			print geometry[geo1]['d'],P1,T1,vapQ1
		counter+=1
	system.simData['vars']=sim_vars
	return system

def forward_steady_calculate(system, tol_P=10,flow='b_b',VLE='ideal',tol_T=1):
	system=forward_steady_itterate(system,tol_P=tol_P,flow=flow,VLE=VLE)
	P_out=system.simData['vars'][-1]['P']
	error=tol_P*2
	it_counter=0
	devs=32
	while error>tol_P:
		print P_out
		system=forward_steady_itterate(system,tol_P=tol_P,flow=flow,VLE=VLE,devs=devs,tol_T=tol_T)
		P_out_new=system.simData['vars'][-1]['P']
		error=abs(P_out_new-P_out)
		P_out=P_out_new
		devs=devs*2
		it_counter+=1
	system.transSimData=[{'t':0,'simData':system.simData}]
	return (system,it_counter)