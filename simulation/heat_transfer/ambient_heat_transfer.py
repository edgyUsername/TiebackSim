import math
from scipy.optimize import brentq
ln=math.log
exp=math.exp
def calculate_T(D,L,massFlow,Cp_l,Cp_g,vapQ,T0,T_out,U):
	A=math.pi*D*L
	Cp_g=0 if Cp_g==None else Cp_g
	Cp_overall=vapQ*Cp_g+(1-vapQ)*Cp_l
	T1=(T0 - T_out)*math.exp(-U*A/(massFlow*Cp_overall))+T_out
	return T1
# def Q_wall(T,D,L,T0,T_out,U):												#-
# 	A = math.pi * D * L
# 	Q=U*A*(T-T0)/ln((T-T_out)/(T0-T_out))
# 	return Q
# def Q_fluid(T,massFlow,Cp_l,Cp_g,vapQ,T0):									#-
# 	Cp_g = 0 if Cp_g == None else Cp_g
# 	Cp_overall = vapQ * Cp_g + (1 - vapQ) * Cp_l
# 	Q=massFlow*Cp_overall*(T-T0)
# 	return Q
# def Q_internal(h1,h0,massFlow):												#-
# 	return (h1-h0)*massFlow

def minimize_q(D,L,massFlow,Cp_l,Cp_g,vapQ,T0,T1,T_out,U,h1,h0):
	Cp= vapQ * Cp_g + (1 - vapQ) * Cp_l
	A = math.pi * D * L
	def fun(x):
		return -x+T_out + (T0 - T_out) * exp(-(x - T0) * U * A / ((h1 - h0) * massFlow))
	tolT=0.05
	error=tolT+1
	print 'TTTTTTTTTTTTTTTTTTTTTTTTTTTT\n',T1
	# T=brentq(fun,max(T0,T1,T_out)+50,min(T0,T1,T_out)-50)
	# T=2*((h1-h0)*massFlow/(T1-T0)+T_out)/(1+2*(h1-h0)/(T1-T0))
	T=T1+((h1-h0)*massFlow+U*A*(T1-T0)/ln((T1-T_out)/(T0-T_out)))/(massFlow*Cp)
	print T1,T
	# while error>tolT:
	# 	T=fun(T1)
	# 	print T
	# 	error=abs(T-T1)
	# 	T1=T
	print 'LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL'
	return T
