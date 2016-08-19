import math
from scipy.optimize import minimize_scalar
def ln(x):
	if abs(x)>1e-15:
		return math.log(x)
	else:
		return -1e17
def get_wilson_lnK(pvt,T,P):
	K={}
	for i in pvt:
		K[i]=ln(pvt[i]['params']['p_c']*1e5/P)+5.373*(1+pvt[i]['params']['acc'])*(1-pvt[i]['params']['t_c']/(T+273.14))
	return K
def rachford_rice_g(vapour_frcn,pvt,K):
	#molar fracn
	g=0
	for i in pvt:

		try:
			z,K_i=pvt[i]['comp'],math.exp(K[i])
		except OverflowError:
			z,K_i=pvt[i]['comp'],1e15
			print 'overflow error'
		try:
			g+=z*(K_i-1)/(1- vapour_frcn*(1- K_i))
		except ZeroDivisionError:
			g+=1e15
			print 'ZeroDivision error'
	return g
# def find_vapor_frcn(pvt,K):
# 	print'K: ', K
# 	g_0=rachford_rice_g(0,pvt,K)
# 	comp_v={}
# 	comp_l={}
# 	V=-1
# 	if g_0<=0:	#bubble point
# 		V=0.0
# 		for i in pvt:
# 			comp_v[i]=pvt[i]['comp']*math.exp(K[i])
# 			comp_l[i]=pvt[i]['comp']
# 	else:
# 		g_1=rachford_rice_g(1,pvt,K)
# 		if g_1>=0: #dew point
# 			V=1.0
# 			for i in pvt:
# 				comp_v[i]=pvt[i]['comp']
# 				comp_l[i]=pvt[i]['comp']/math.exp(K[i])
# 	if V==-1:
# 		tol=0.00001
# 		e=tol+1
# 		Va=0.0
# 		ga=g_0
# 		Vb=1.0
# 		gb=g_1
# 		itc=0
# 		while e>tol:
# 			Vc=(Va+Vb)/2
# 			gc=rachford_rice_g(Vc,pvt,K)
# 			if gc*ga>0:
# 				Va=Vc
# 				ga=gc
# 			else:
# 				Vb=Vc
# 				gb=gc
# 			e=abs(Va-Vb)
# 			itc+=1
# 		V=Vc
# 		for i in pvt:
# 			comp_v[i]=pvt[i]['comp']*math.exp(K[i])/(1-V*(1- math.exp(K[i])))
# 			comp_l[i]=pvt[i]['comp']/(1-V*(1- math.exp(K[i])))
# 	return (V,comp_v,comp_l)
def find_vapor_frcn(pvt,K):
	V=minimize_scalar(lambda x:rachford_rice_g(x,pvt,K)**2,bracket=[0,1],bounds=[0,1],method='bounded').x
	comp_l={}
	comp_v={}
	for i in pvt:
		comp_v[i] = pvt[i]['comp'] * math.exp(K[i]) / (1 - V * (1 - math.exp(K[i])))
		comp_l[i] = pvt[i]['comp'] / (1 - V * (1 - math.exp(K[i])))
	return (V,comp_v,comp_l)
def cubic_solver(a,b,c,d):
	"""returns real roots"""
	roots=[]
	delta=18*a*b*c*d-4*b*b*b*d+b*b*c*c-4*a*c*c*c-27*a*a*d*d
	delta_0=b*b-3*a*c
	delta_1=2*b*b*b-9*a*b*c+27*a*a*d
	if delta==0:
		if delta_0==0:
			roots.append(-b/(3*a))
		else:
			roots.append((9*a*d-b*c)/(2*delta_0))
			roots.append((4*a*b*c-9*a*a*d-b*b*b)/(a*delta_0))
	else:
		s=delta_1*delta_1-4*delta_0*delta_0*delta_0
		s_comp=s<0
		if s_comp:
			S=(delta_1-math.sqrt(abs(s))*1j)/2.0
			S=S if not S.real==0 else (delta_1+math.sqrt(abs(s))*1j)/2.0
		else:
			S=(delta_1-math.sqrt(s))/2.0
			S= S if not S==0 else (delta_1+math.sqrt(s))/2.0
			S=S+0j
		C=[]
		r=math.sqrt(S.real*S.real+S.imag*S.imag)
		if S.real==0:
			if S.imag>=0:
				theta=math.pi/2
			else:
				theta=-math.pi/2
		else:
			theta=math.atan(S.imag/S.real)
		if S.real<0:
			if S.imag<0:
				arg=-math.pi+theta
			else:
				arg=math.pi+theta
		else:
			if S.imag<0:
				arg=theta
			else:
				arg=theta
		r_root=math.pow(r,1.0/3)
		C.append(r_root*(math.cos(arg/3)+1j*math.sin(arg/3)))
		C.append(r_root*(math.cos(arg/3+math.pi*2/3)+1j*math.sin(arg/3+math.pi*2/3)))
		C.append(r_root*(math.cos(arg/3+math.pi*4/3)+1j*math.sin(arg/3+math.pi*4/3)))
		reals=0
		tol=1e-8
		root=[]
		reals=0
		for i in C:
			x_k=-(b+i+delta_0/i)/(3*a)
			error=abs(x_k.imag)/((x_k.real**2+x_k.imag**2)**.5)
			if x_k.imag==0 or error<tol:
				reals+=1
				root.append(x_k.real)
		roots=root
	roots.sort()
	return roots
def van_der_waal_coefs(pvt,T,P,comp=None):
	"""returns dict of a and b params for mixtures and comps"""
	if not comp:
		comp={}
		for i in pvt:
			comp[i]=pvt[i]['comp']
	#a_m,b_m,comp_params:{name:{a_i,b_i,alpha_i,m_i,a_ci}}
	comp_params={}
	for i in comp:
		pvt_i=pvt[i]['params']
		if pvt_i['acc']<=.49:
			m_i=.37464+1.54226*pvt_i['acc']-.26992*(pvt_i['acc']**2)
		else:
			m_i=.379642*(1.48503-(.164423-1.016666*pvt_i['acc'])*pvt_i['acc'])*pvt_i['acc']
		alpha_i=(1+m_i*(1- math.sqrt((T+273.14)/pvt_i['t_c'])))**2
		a_ci=.457235*((8.314*pvt_i['t_c'])**2)/(pvt_i['p_c']*1e5)
		a_i=a_ci*alpha_i
		b_i=.077796*8.314*pvt_i['t_c']/(1e5*pvt_i['p_c'])
		comp_params[i]={
		'm_i':m_i,
		'alpha_i':alpha_i,
		'a_ci':a_ci,
		'a_i':a_i,
		'b_i':b_i
		}
	a_m=0
	b_m=0
	for i in comp:
		b_i=comp_params[i]['b_i']
		b_m+=comp[i]*b_i
		a_ij=0
		b_ij=0
		a_i=comp_params[i]['a_i']
		x_i=comp[i]
		for j in comp:
			x_j=comp[j]
			a_j=comp_params[j]['a_i']
			b_j=comp_params[j]['b_i']
			a_ij+=x_j*math.sqrt(a_i*a_j)*(1-pvt[i]['k'][j])
			b_ij+=x_j*(b_j+b_i)/2.
		a_m+=a_ij*x_i
		comp_params['a_ij']=a_ij
		comp_params['b_ij']=b_ij
	return {'a_m':a_m,'b_m':b_m,'comp_params':comp_params}
def solve_PR_for_Z(pvt,T,P,comp=None):
	""" if comp is provided, will give phase specific V.
		returns list of V roots,EoS params
	"""
	R=8.314
	if not comp:
		comp={}
		for i in pvt:
			comp[i]=pvt[i]['comp']
	eos_params=van_der_waal_coefs(pvt,T,P,comp=comp)
	a=eos_params['a_m']
	b=eos_params['b_m']
	T,R=T+273.15,8.314
	A_1,B_1=a*P/((R*T)**2),b*P/(R*T)
	A,B,C,D=1,(B_1-1),(A_1- 2*B_1 -3*(B_1**2)),-(A_1*B_1-(B_1**2)-(B_1**3))
	Z= cubic_solver(A,B,C,D)
	V=[R*T*i/P for i in Z]
	# if len(V)==0:
	# 	print Tp,P,comp
	return (Z,eos_params)
def fug_minimum_gibbs(pvt,T,P,n,comp=None,returnV=False,phase="None"):
	"""calculates fug coefs that correspond to minimum gibbs energy"""
	if not comp:
		comp={}
		for i in pvt:
			comp[i]=pvt[i]['comp']
	Z_roots,eos_params=solve_PR_for_Z(pvt,T,P,comp=comp)
	a,b=eos_params['a_m'],eos_params['b_m']
	R,T=8.314,T+273.15
	Z_l=Z_roots[-1]		#gas like phase
	Z_h=Z_roots[0]			#liquid like phase
	V_l0=Z_l*R*T/P
	V_h0=Z_h*R*T/P
	V_l,V_h=V_l0 if V_l0>b else (1+1e-4)*b,V_h0 if V_h0>b else (1+1e-4)*b
	Z_l,Z_h=Z_l if V_l0>b else V_l*P/(R*T),Z_h if V_h0>b else V_h*P/(R*T)
	lnfug_l,lnfug_h={},{}
	A,B=a*P/((R*T)**2),b*P/(R*T)
	for i in pvt:
		if a==0 or b==0:
			if phase=='light':
				ln_fug_l =100
				ln_fug_h=-100
			elif phase=='heavy':
				ln_fug_l =-100
				ln_fug_h =100
		else:
			AA=2*sum([comp[j]*math.sqrt(eos_params['comp_params'][j]['a_i']*eos_params['comp_params'][i]['a_i']) for j in pvt])/a
			BB=eos_params['comp_params'][i]['b_i']/b
			ln_fug_l=BB*(Z_l-1)-ln(Z_l-B)-A*(AA-BB)*ln((Z_l+((2**.5)+1)*B)/((Z_l-((2**.5)-1)*B)))/((2**1.5)*B)
			ln_fug_h=BB*(Z_h-1)-ln(Z_h-B)-A*(AA-BB)*ln((Z_h+((2**.5)+1)*B)/((Z_h-((2**.5)-1)*B)))/((2**1.5)*B)
		lnfug_l[i]=ln_fug_l
		lnfug_h[i]=ln_fug_h

	if phase=='light':
		if returnV:
			return lnfug_l,V_l
		else:
			return lnfug_l
	elif phase=='heavy':
		if returnV:
			return lnfug_h,V_h
		else:
			return lnfug_h
	else:
		fug_l=0
		fug_h=0
		for i in comp:
			fug_l+=ln(comp[i])+lnfug_l[i]
			fug_h+=ln(comp[i])+lnfug_h[i]
		if fug_l<fug_h:
			fug,V=lnfug_l,V_l
		else:
			fug,V=lnfug_h,V_h
		if returnV:
			return fug,V
		else:
			return fug
