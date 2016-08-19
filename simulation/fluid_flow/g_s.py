import math
log=math.log10
# def calcDWf(roughness,ID,Re):
# 	"""
# 	goudar-sonar eq for D-W f
# 	"""
# 	a=2/math.log(10)
# 	b=roughness/(ID*3.7)
# 	d=math.log(10)*Re/(5.02)
# 	s=b*d+math.log(d)
# 	q=s**(s/(s+1))
# 	g=b*d+math.log(d/q)
# 	z=math.log(q/g)
# 	dLA=z*g/(g+1)
# 	dCFA=dLA*(1+(z/2)/(((g+1)**2)+(z/3.0)*(2*g-1)))
# 	f=(a*(math.log(d/q)+dCFA))**-2
# 	return f
def calcDWf(roughness,ID,Re):
	if Re<2000:
		return 64./Re
	else:
		F=lambda F: -2*log(roughness/(3.7*ID)+2.51*F/Re)
		F0=1/math.sqrt(.0001)
		tol=1/math.sqrt(1e-7)
		error=tol+1
		while error>tol:
			F1=F(F0)
			error=abs(F1-F0)
			F0=F1
		f=1/(F1**2)
		return f
def calcFf(roughness,ID,Re):
	return calcDWf(roughness,ID,Re)/4