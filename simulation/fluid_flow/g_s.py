import math

def calcDWf(roughness,ID,Re):
	"""
	goudar-sonar eq for D-W f
	"""
	print roughness,ID,Re
	a=2/math.log(10)
	b=roughness/(ID*3.7)
	d=math.log(10)*Re/(5.02)
	s=b*d+math.log(d)
	print b,d,s
	q=s**(s/(s+1))
	g=b*d+math.log(d/q)
	z=math.log(q/g)
	dLA=z*g/(g+1)
	dCFA=dLA*(1+(z/2)/(((g+1)**2)+(z/3.0)*(2*g-1)))
	f=(a*(math.log(d/q)+dCFA))**-2
	return f

def calcFf(roughness,ID,Re):
	return calcDWf(roughness,ID,Re)/4