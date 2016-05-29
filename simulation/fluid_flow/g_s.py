import math

def calcDWf(roughness,ID,Re):
	"""
	goudar-sonar eq for D-W f
	"""
	a=2/math.log(10,math.e)
	b=roughness/(ID*3.7)
	d=math.log(10,math.e)*Re/(5.02)
	s=b*d+math.log(d,math.e)
	q=math.pow(s,s/(s+1))
	g=b*d+math.log(d/q,math.e)
	z=math.log(q/g,math.e)
	dLA=z*g/(g+1)
	dCFA=dLA*(1+(z/2)/(math.pow((g+1),2)+(z/3.0)*(2*g-1)))
	f=math.pow(a*(math.log(d/q,math.e)+dCFA),-2)
	return f