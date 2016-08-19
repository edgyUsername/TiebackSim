import sys
import math
from os import path
sys.path.append( path.dirname( path.dirname(path.dirname( path.abspath(__file__)) ) ) )
#import settings
from simulation.fluid_flow import g_s

def ln(x):
	return math.log(x,math.e)

bb_params={
	'segregated':{
		'a':0.98,
		'b':0.4846,
		'c':0.0868,
		'd':0.011,
		'e':-3.768,
		'f':3.539,
		'g':-1.614
	},
	'intermittent':{
		'a':0.845,
		'b':0.5351,
		'c':0.0173,
		'd':2.96,
		'e':0.305,
		'f':-0.4473,
		'g':0.0978
	},
	'distributed':{
		'a':1.065,
		'b':0.5351,
		'c':0.0609,
	},
	'downhill':{
		'd':4.7,
		'e':-0.3692,
		'f':0.1244,
		'g':-0.5056
	},
}


def estimateSigma (r_l,T,P):
	T_F=T*1.8+32
	if T_F<68:
		T_F=68
	elif T_F>100:
		T_F=100
	print r_l
	sg=r_l/999.016
	API=141.5/sg-131.5
	print '############# API ##############'
	print '              ',API
	sigma68=39-0.2571*API
	sigma100=37.5-0.2571*API
	print '              ',sigma68,sigma100
	sigma=sigma68-(T_F-68)*(sigma68-sigma100)/32
	print'              ',sigma
	print '################################'
	p=P/6894.757
	sigma=sigma*math.exp(-8.6306e-4*p)
	return sigma/1000

def calcPressureDrop(P0,m,vapQ,r_g,r_l,m_g,m_l,z,roughness,ID,L,T):
	Q_l=m*(1-vapQ)/r_l
	Q_g=m*vapQ/r_g
	H=Q_l/(Q_g+Q_l)
	A=math.pi*math.pow(ID,2)/4
	v_sg=Q_g/A
	v_sl=Q_l/A
	v_m=v_sg+v_sl
	nFR=math.pow(v_m,2)/(9.81*ID)
	#lambda_l=v_l/v_m
	lambda_l=H

	L1=316*math.pow(lambda_l,0.302)
	L2=0.0009252*math.pow(lambda_l,-2.4684)
	L3=0.1*math.pow(lambda_l,-1.4516)
	L4=0.5*math.pow(lambda_l,-6.738)

	regime=None
	if (lambda_l<0.01 and nFR<L1) or (lambda_l>=0.01 and nFR<L2):
		regime='segregated'
	elif lambda_l>=0.01 and L2<nFR<=L3:
		regime='transition'
	elif (0.01<=lambda_l<4 and L3<nFR<=L1) or (lambda_l>=0.4 and L3<nFR<=L4):
		regime='intermittent'
	elif (lambda_l<0.4 and nFR>=L1) or (lambda_l>=0.4 and nFR>L4):
		regime='distributed'
	else:
		assert False
	theta=math.asin(max(min(z/L,1),-1))
	sigma=estimateSigma(r_l,T,P0)
	nVl=v_sl*math.pow(r_l/(9.81*sigma*0.001),0.25)
	if regime=='segregated':
		a=bb_params['segregated']['a']
		b=bb_params['segregated']['b']
		c=bb_params['segregated']['c']
		yl0=a*math.pow(lambda_l,b)/math.pow(nFR,c)
		yl0=yl0 if yl0>=lambda_l else lambda_l
		if z<=0:
			d=bb_params['downhill']['d']
			e=bb_params['downhill']['e']
			f=bb_params['downhill']['f']
			g=bb_params['downhill']['g']
		elif z>0:
			d=bb_params['segregated']['d']
			e=bb_params['segregated']['e']
			f=bb_params['segregated']['f']
			g=bb_params['segregated']['g']
		else:
			assert False
		C=(1-lambda_l)*math.log(d*math.pow(lambda_l,e)*math.pow(nVl,f)*math.pow(nFR,g),math.e)
		psi=1+C*(math.sin(1.8*theta)-0.333*math.pow(math.sin(1.8*theta),3))
		yl=yl0*psi
	elif regime=='intermittent':
		a=bb_params['intermittent']['a']
		b=bb_params['intermittent']['b']
		c=bb_params['intermittent']['c']
		yl0=a*math.pow(lambda_l,b)/math.pow(nFR,c)
		yl0=yl0 if yl0>=lambda_l else lambda_l
		if z<=0:
			d=bb_params['downhill']['d']
			e=bb_params['downhill']['e']
			f=bb_params['downhill']['f']
			g=bb_params['downhill']['g']
		elif z>0:
			d=bb_params['intermittent']['d']
			e=bb_params['intermittent']['e']
			f=bb_params['intermittent']['f']
			g=bb_params['intermittent']['g']
		else:
			assert False
		C=(1-lambda_l)*math.log(d*math.pow(lambda_l,e)*math.pow(nVl,f)*math.pow(nFR,g),math.e)
		psi=1+C*(math.sin(1.8*theta)-0.333*math.pow(math.sin(1.8*theta),3))
		yl=yl0*psi
	elif regime=='distributed':
		a=bb_params['distributed']['a']
		b=bb_params['distributed']['b']
		c=bb_params['distributed']['c']
		yl0=a*math.pow(lambda_l,b)/math.pow(nFR,c)
		yl0=yl0 if yl0>=lambda_l else lambda_l
		if z<=0:
			d=bb_params['downhill']['d']
			e=bb_params['downhill']['e']
			f=bb_params['downhill']['f']
			g=bb_params['downhill']['g']
			C=(1-lambda_l)*math.log(d*math.pow(lambda_l,e)*math.pow(nVl,f)*math.pow(nFR,g),math.e)
			psi=1+C*(math.sin(1.8*theta)-0.333*math.pow(math.sin(1.8*theta),3))
			yl=yl0*psi
		elif z>0:
			yl=yl0
		else:
			assert False

	elif regime=='transition':
		A=(L3-nFR)/(L3-L2)
		B=1-A
		####segregated####
		a=bb_params['segregated']['a']
		b=bb_params['segregated']['b']
		c=bb_params['segregated']['c']
		yls0=a*math.pow(lambda_l,b)/math.pow(nFR,c)
		yls0=yls0 if yls0>=lambda_l else lambda_l
		if z<=0:
			d=bb_params['downhill']['d']
			e=bb_params['downhill']['e']
			f=bb_params['downhill']['f']
			g=bb_params['downhill']['g']
		elif z>0:
			d=bb_params['segregated']['d']
			e=bb_params['segregated']['e']
			f=bb_params['segregated']['f']
			g=bb_params['segregated']['g']
		else:
			assert False
		C=(1-lambda_l)*math.log(d*math.pow(lambda_l,e)*math.pow(nVl,f)*math.pow(nFR,g),math.e)
		psi=1+C*(math.sin(1.8*theta)-0.333*math.pow(math.sin(1.8*theta),3))
		yls=yls0*psi
		###intermitent###
		a=bb_params['intermittent']['a']
		b=bb_params['intermittent']['b']
		c=bb_params['intermittent']['c']
		yli0=a*math.pow(lambda_l,b)/math.pow(nFR,c)
		yli0=yli0 if yli0>=lambda_l else lambda_l
		if z<=0:
			d=bb_params['downhill']['d']
			e=bb_params['downhill']['e']
			f=bb_params['downhill']['f']
			g=bb_params['downhill']['g']
		elif z>0:
			d=bb_params['intermittent']['d']
			e=bb_params['intermittent']['e']
			f=bb_params['intermittent']['f']
			g=bb_params['intermittent']['g']
		else:
			assert False
		C=(1-lambda_l)*math.log(d*math.pow(lambda_l,e)*math.pow(nVl,f)*math.pow(nFR,g),math.e)
		psi=1+C*(math.sin(1.8*theta)-0.333*math.pow(math.sin(1.8*theta),3))
		yli=yli0*psi
		####combine###
		yl=A*yls+B*yli
	yl=1.0 if yl>1.0 else yl
	### gravitational potential
	r_m_is=yl*r_l+(1-yl)*r_g
	dp_PE=z*r_m_is*9.81
	### frictional
	r_m=r_l*lambda_l+r_g*(1-lambda_l)
	m_m=m_l*lambda_l+m_g*(1-lambda_l)
	Re=r_m*v_m*ID/m_m
	f=g_s.calcDWf(roughness,ID,Re)
	x=lambda_l/math.pow(yl,2)

	if 1<x<1.2:
		S=ln(2.2*x-1.2)
	else:
		S=ln(x)/(-0.0523+3.182*ln(x)-0.8725*math.pow(ln(x),2)+0.01853*math.pow(ln(x),4))
	f_tp=f*math.exp(S)
	dp_f=0.5*f_tp*math.pow(v_m,2)*r_m*L/ID
	####Ek####
	# Ek=v_m*v_sg*r_m/P0
	####Total pressure drop####
	# dP_total=(dp_f+dp_PE)/(1-Ek)
	dP_total = (dp_f + dp_PE)
	return (P0-dP_total,yl,regime)






