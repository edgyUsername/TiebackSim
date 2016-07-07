import math
from simulation.fluid_flow import g_s
def calcPressureDrop(P0,m,r,mu,z,roughness,ID,L,T):
	v=m/(r*ID*ID*math.pi/4)
	Re=4*m/(ID*math.pi*mu)
	f=g_s.calcDWf(roughness,ID,Re)
	P=P0-r*9.81*(z+0.5*f*(v**2)*L/(9.81*ID))
	return P