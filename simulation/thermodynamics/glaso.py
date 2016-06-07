import math
def updateLiquidVisc(m_l0,r_l,T1,T0):
	T0=(T0+273.14)*1.8
	T1=(T1+273.14)*1.8
	API=141.5*1000/r_l-131.5
	a1=10.313*math.log10(T1-460)-36.447
	a0=10.313*math.log10(T0-460)-36.447
	m_l1=m_l0*math.pow((T1-460)/(T0-460),-3.444)*math.pow(math.log10(API),a1-a0)
	return m_l1