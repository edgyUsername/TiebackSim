import math

def calculate_T(D,L,massFlow,Cp_l,Cp_g,vapQ,T0,T_out,U):
	A=math.pi*D*L
	Cp_g=0 if Cp_g==None else Cp_g
	Cp_overall=vapQ*Cp_g+(1-vapQ)*Cp_l
	T1=(T0 - T_out)*math.exp(-U*A/(massFlow*Cp_overall))+T_out
	return T1