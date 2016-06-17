import math

def normalise_comp(comp):
	total=0
	original={}

	for i in comp:
		original[comp[i]]=i
		total+=comp[i] if comp[i]>0 else 1e-7
	for i in comp:
		comp[i],'/',total
		comp[i]=comp[i]/float(total) if comp[i]>0 else 1e-7
	total=0
	for i in comp:
		total+=comp[i]
	if total==0:
		comp[original[max(original.keys())]]=1
	return comp
def ln(x):
	if abs(x)>1e-15:
		return math.log(x)
	else:
		return -1e15
	
def round_sig(x, sig=14):
	return round(x, sig-int(math.floor(math.log10(abs(x))))-1)
def get_wilson_K(pvt,T,P):
	K={}
	for i in pvt:
		K[i]=math.exp(ln(pvt[i]['params']['p_c']*1e5/P)+5.373*(1+pvt[i]['params']['acc'])*(1-pvt[i]['params']['t_c']/(T+273.14)))
	return K
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
		try:
			g+=z*(K_i-1)/(1- vapour_frcn*(1- K_i))
		except ZeroDivisionError:
			g+=1e15
	return g
def find_vapor_frcn(pvt,K):
	g_0=rachford_rice_g(0,pvt,K)
	comp_v={}
	comp_l={}
	V=-1
	if g_0<=0:	#bubble point
		V=0.0
		for i in pvt:
			comp_v[i]=pvt[i]['comp']*math.exp(K[i])
			comp_l[i]=pvt[i]['comp']
	else:
		g_1=rachford_rice_g(1,pvt,K)
		if g_1>=0: #dew point
			V=1.0
			for i in pvt:
				comp_v[i]=pvt[i]['comp']
				comp_l[i]=pvt[i]['comp']/math.exp(K[i])
	if V==-1:
		tol=0.00001
		e=tol+1
		Va=0.0
		ga=g_0
		Vb=1.0
		gb=g_1
		itc=0
		while e>tol:
			Vc=(Va+Vb)/2
			gc=rachford_rice_g(Vc,pvt,K)
			if gc*ga>0:
				Va=Vc
				ga=gc
			else:
				Vb=Vc
				gb=gc
			e=abs((Va-Vb)/Vc)
			itc+=1
		V=Vc
		for i in pvt:
			comp_v[i]=pvt[i]['comp']*math.exp(K[i])/(1-V*(1- math.exp(K[i])))
			comp_l[i]=pvt[i]['comp']/(1-V*(1- math.exp(K[i])))
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
			a_ij+=x_j*math.sqrt(a_i*a_j)
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
	T,R=T+273.16,8.314
	A_1,B_1=a*P/((R*T)**2),b*P/(R*T)
	A,B,C,D=1,(B_1-1),(A_1- 2*B_1 -3*(B_1**2)),-(A_1*B_1-(B_1**2)-(B_1**3))
	Z= cubic_solver(A,B,C,D)
	V=[R*T*i/P for i in Z]
	if len(V)==0:
		print T,P,comp
	return (Z,eos_params)
def ln_fug_coef(pvt,T,P,comp=None,phase="light",returnV=False):
	"""returns dict of ln of fugacity coefficients"""
	if not comp:
		comp={}
		for i in pvt:
			comp[i]=pvt[i]['comp']
	V_roots,eos_params=solve_PR_for_V(pvt,T,P,comp=comp)
	if phase=='light':
		V_m=V_roots[-1]
	elif phase=='heavy':
		V_m=V_roots[0]
	else:
		assert False,"phase must be either 'heavy' or 'light'"
	T_k=T+273.14
	V_ideal=8.314*T_k/P
	ln_fug_coef={}
	a_m,b_m=eos_params['a_m'],eos_params['b_m']
	V_m = b_m*1.0001 if V_m<=b_m else V_m
	for i in comp:
		i_params=eos_params['comp_params'][i]
		a_i,b_i=i_params['a_i'],i_params['b_i']
		ln_fug_i=(P*V_m/(8.314*T_k)-1
			+ln(V_ideal/(V_m- b_m))
			+(b_i*a_m-2*math.sqrt(a_i*a_m)*b_m)*ln((V_m+b_m+math.sqrt(2)*b_m)/(V_m+b_m-math.sqrt(2)*b_m))/(2*math.sqrt(2)*8.314*T_k*b_m*b_m)
			+(b_i- b_m)/(V_m- b_m)
			+(b_i- b_m)*a_m*V_m/(8.314*T_k*(2*b_m*b_m-math.pow(V_m+b_m,2))))
		ln_fug_coef[i]=ln_fug_i
	return ln_fug_coef if not returnV else (ln_fug_coef,V_m)
def get_eos_lnK(pvt,T,P,comp_l,comp_h,V):
	lnK={}
	ln_fug_coef_l=fug_minimum_gibbs(pvt,T,P,V,comp=comp_l,phase='light')
	ln_fug_coef_h=fug_minimum_gibbs(pvt,T,P,(1-V),comp=comp_h,phase='heavy')
	
	for i in pvt:
		lnK[i]=ln_fug_coef_h[i]-ln_fug_coef_l[i]
	return lnK
def successive_substitution(pvt,T,P,comp_l,comp_h,V):
	"""returns V, new comp_l,new comp_h"""
	lnK=get_eos_lnK(pvt,T,P,comp_l,comp_h,V)
	print 'successive sub: ',lnK
	return find_vapor_frcn(pvt,lnK)
def tpd_analysis(pvt,T,P,comp,comp_l,comp_h,V,output='stability'):
	
	ln_fug_spec=fug_minimum_gibbs(pvt,T,P,1,comp=comp)
	ln_fug_h=fug_minimum_gibbs(pvt,T,P,(1-V),comp=comp_h,phase="heavy")
	ln_fug_l=fug_minimum_gibbs(pvt,T,P,V,comp=comp_l,phase="light")
	tpd_l=0
	tpd_h=0
	for i in comp:
		x_h=comp_h[i]
		x_l=comp_l[i]
		x=comp[i]
		tpd_l+=x_l*(ln(x_l)+ln_fug_l[i]-ln_fug_spec[i]-ln(x)) if x_l>0 else 0
		tpd_h+=x_h*(ln(x_h)+ln_fug_h[i]-ln_fug_spec[i]-ln(x)) if x_h>0 else 0	
	delta_G_nrt=(1-V)*tpd_h+V*tpd_l
	stable= (delta_G_nrt>=0) and (tpd_l>=0) and (tpd_h>=0)
	if output=='stability':
		return stable
	elif output=='tpd':
		return tpd_l,tpd_h,delta_G_nrt
	else:
		assert False, "output must be either tpd or stability"
def fug_minimum_gibbs(pvt,T,P,n,comp=None,returnV=False,phase="None"):
	"""calculates fug coefs that correspond to minimum gibbs energy"""
	if not comp:
		comp={}
		for i in pvt:
			comp[i]=pvt[i]['comp']
	Z_roots,eos_params=solve_PR_for_Z(pvt,T,P,comp=comp)
	a,b=eos_params['a_m'],eos_params['b_m']
	#D,B=a*n*n,b*n
	R,T=8.314,T+273.14
	Z_l=Z_roots[-1]		#gas like phase
	Z_h=Z_roots[0]			#liquid like phase
	V_l0=Z_l*R*T/P
	V_h0=Z_h*R*T/P
	V_l,V_h=V_l0 if V_l0>b else (1+1e-4)*b,V_h0 if V_h0>b else (1+1e-4)*b
	Z_l,Z_h=Z_l if V_l0>b else V_l*P/(R*T),Z_h if V_h0>b else V_h*P/(R*T)
	lnfug_l,lnfug_h={},{}
	A,B=a*P/((R*T)**2),b*P/(R*T)
	for i in pvt:
		#### derivative ###
		# Bi=eos_params['comp_params'][i]['b_i']
		# sum_a=0
		# for j in pvt:
		# 	sum_a+=comp[j]*math.sqrt(eos_params['comp_params'][j]['a_i']*eos_params['comp_params'][i]['a_i'])
		# Di=2*n*sum_a
		# ### light ###
		# Fn_l=ln(V_l)-ln(V_l-B)-D*B*B*n*ln(1+B*n/V_l-(B*n/V_l)*(B*n/V_l))/(R*T*V_l*V_l*(B-B*B*n*n/(V_l*V_l))*(B-B*B*n*n/(V_l*V_l)))-D*(2*B/V_l-2*B*B*n/(V_l*V_l))/(2*R*T*(B-B*B*n*n/(V_l*V_l))*(1+2*B*n/V_l-B*B*n*n/(V_l*V_l)))
		# FB_l=D*(1-2*n*n*B/(V_l*V_l))*ln(1+2*n*B/V_l-n*n*B*B/(V_l*V_l))/(2*R*T*(B-n*n*B*B/(V_l*V_l))*(B-n*n*B*B/(V_l*V_l)))-D*(2*n/V_l-2*n*n*B/(V_l*V_l))/(2*R*T*(B-n*n*B*B/(V_l*V_l))*(1+2*n*B/V_l-n*n*B*B/(V_l*V_l)))+n/(V_l-B)
		# FD_l=-ln(1+2*B*n/V_l-B*B*n*n/(V_l*V_l))/(2*R*T*(B-B*B*n*n/(V_l*V_l)))
		# ln_fug_l=Fn_l+FB_l*Bi+FD_l*Di-ln(Z_l)
		# lnfug_l[i]=ln_fug_l
		# ### heavy ###
		# Fn_h=ln(V_h)-ln(V_h-B)-D*B*B*n*ln(1+B*n/V_h-(B*n/V_h)*(B*n/V_h))/(R*T*V_h*V_h*(B-B*B*n*n/(V_h*V_h))*(B-B*B*n*n/(V_h*V_h)))-D*(2*B/V_h-2*B*B*n/(V_h*V_h))/(2*R*T*(B-B*B*n*n/(V_h*V_h))*(1+2*B*n/V_h-B*B*n*n/(V_h*V_h)))
		# FB_h=D*(1-2*n*n*B/(V_h*V_h))*ln(1+2*n*B/V_h-n*n*B*B/(V_h*V_h))/(2*R*T*(B-n*n*B*B/(V_h*V_h))*(B-n*n*B*B/(V_h*V_h)))-D*(2*n/V_h-2*n*n*B/(V_h*V_h))/(2*R*T*(B-n*n*B*B/(V_h*V_h))*(1+2*n*B/V_h-n*n*B*B/(V_h*V_h)))+n/(V_h-B)
		# FD_h=-ln(1+2*B*n/V_h-B*B*n*n/(V_h*V_h))/(2*R*T*(B-B*B*n*n/(V_h*V_h)))
		# ln_fug_h=Fn_h+FB_h*Bi+FD_h*Di-ln(Z_h)
		# lnfug_h[i]=ln_fug_h
		#### method2 ####
		if Z_l<=0 or Z_h<=0:
			print Z_roots
		AA=2*sum([math.sqrt(eos_params['comp_params'][j]['a_i']*eos_params['comp_params'][i]['a_i']) for j in pvt])/a
		BB=eos_params['comp_params'][i]['b_i']/b
		ln_fug_l=BB*(Z_l-1)-ln(Z_l-B)-A*(AA-BB)*ln((Z_l+((2**.5)+1)*B)/((Z_l-((2**.5)-1)*B)))/((2**1.5)*B)
		ln_fug_h=BB*(Z_h-1)-ln(Z_h-B)-A*(AA-BB)*ln((Z_h+((2**.5)+1)*B)/((Z_h-((2**.5)-1)*B)))/((2**1.5)*B)
		lnfug_l[i]=ln_fug_l
		lnfug_h[i]=ln_fug_h
		####### TEST TEST TEST ####
		# e=1e-5
		# new_comp1={}
		# new_comp2={}
		# for z in comp:
		# 	new_comp2[z]=comp[z] if not z==i else comp[z]+e
		# 	new_comp1[z]=comp[z] if not z==i else comp[z]-e
		# (V_1,eos1),(V_2,eos2)=solve_PR_for_V(pvt,(T-273.14),P,comp=normalise_comp(new_comp1)),solve_PR_for_V(pvt,(T-273.14),P,comp=normalise_comp(new_comp2))
		# n1,n2=1-e,1+e
		# V_1,V_2=V_1[0]*n1,V_2[0]*n2
		# dF=Fn_h+FB_h*Bi+FD_h*Di
		# delta_F=(n2*ln(V_2)-ln(V_2-eos2['b_m']*n2)-n2*n2*eos2['a_m']*(ln(V_2+eos2['b_m']*n2)-ln(V_2))/(eos2['b_m']*n2*R*T)-(n1*ln(V_1)-ln(V_1-eos1['b_m']*n1)-n1*n1*eos1['a_m']*(ln(V_1+eos1['b_m']*n1)-ln(V_1))/(eos1['b_m']*n1*R*T)))/(2*e)
	# lnfug_l,V_l=ln_fug_coef(pvt,T,P,comp=comp,phase="light",returnV=True)
	# lnfug_h,V_h=ln_fug_coef(pvt,T,P,comp=comp,phase="heavy",returnV=True)
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

def check_phase_stability(pvt,T,P,comp=None):
	"""returns stability(bool),V,comp_l,comp_h"""
	lnK_wilson=get_wilson_lnK(pvt,T,P)
	if not comp:
		comp={}
		for i in pvt:
			comp[i]=pvt[i]['comp']
	print lnK_wilson
	########### estimate VLE and stab ########
	V,comp_v,comp_l=find_vapor_frcn(pvt,lnK_wilson)
	it=3
	it_count=0
	stable=False
	comp_v,comp_l=normalise_comp(comp_v),normalise_comp(comp_l)
	########### successive substitution #######
	while not stable and it_count<it:
		V,comp_v,comp_l=successive_substitution(pvt,T,P,comp_v,comp_l,V)
		stable= V==0 or V==1
		it_count+=1
	comp_v,comp_l=normalise_comp(comp_v) if V==0 else comp_v,normalise_comp(comp_l) if V==1 else comp_l
	# works up to here
	########### tangent plane analysis ########
	stable=tpd_analysis(pvt,T,P,comp,comp_v,comp_l,V)
	########### detailed stability analysis ########
	# if stable:
	# 	V_trial,comp_v_min,comp_l_min=find_vapor_frcn(pvt,lnK_wilson)
	# 	comp_v_trial,comp_l_trial=comp_v_min,comp_l_min
	# 	ln_fug_v_spec=fug_minimum_gibbs(pvt,T,P,comp=comp)
	# 	ln_fug_l_spec=fug_minimum_gibbs(pvt,T,P,comp=comp)
	# 	tpd_min_v,tpd_min_l,delta_G_min=tpd_analysis(pvt,T,P,comp,normalise_comp(comp_v_min),normalise_comp(comp_l_min),V_trial,output='tpd')
	# 	it,c,c2=10,0,0
	# 	second_counter=False
	# 	while c<10 and c2<3:
	# 		try:
	# 			ln_fug_v_trial=fug_minimum_gibbs(pvt,T,P,comp=comp_v_trial)
	# 			ln_fug_l_trial=fug_minimum_gibbs(pvt,T,P,comp=comp_l_trial)
	# 			stop=True
	# 			stop2=True
	# 			for i in comp:
	# 				try:
	# 					fug_i_v_spec=math.exp(ln_fug_v_spec[i])
	# 				except OverflowError:
	# 					fug_i_v_spec=1e80
	# 				try:
	# 					fug_i_l_spec=math.exp(ln_fug_l_spec[i])
	# 				except OverflowError:
	# 					fug_i_l_spec=1e80
	# 				try:
	# 					fug_i_v=math.exp(ln_fug_v_trial[i])
	# 				except OverflowError:
	# 					fug_i_v=1e80
	# 				try:
	# 					fug_i_l=math.exp(ln_fug_l_trial[i])
	# 				except OverflowError:
	# 					fug_i_l=1e80
	# 				x_i_spec=comp[i]
	# 				x_v_old=comp_v_trial[i]
	# 				x_l_old=comp_l_trial[i]
	# 				comp_v_trial[i]=x_i_spec*fug_i_v_spec/fug_i_v if fug_i_v!=0 else 1e80
	# 				comp_l_trial[i]=x_i_spec*fug_i_l_spec/fug_i_l if fug_i_l!=0 else 1e80
	# 				if abs(comp_v_trial[i]- x_v_old)>1e-5 or abs(comp_l_trial[i]- x_l_old)>1e-5:
	# 					stop=False
	# 				else:
	# 					stable=False
	# 				if abs(comp_v_trial[i]- x_i_spec)>1e-5 or abs(comp_l_trial[i]- x_i_spec)>1e-5:
	# 					stop2=False
	# 		except:
	# 			pass
	# 		if stop or stop2:
	# 			break
	# 		if not second_counter:
	# 			tpd_v_trial,tpd_l_trial,delta_G_trial=tpd_analysis(pvt,T,P,comp,normalise_comp(comp_v_trial),normalise_comp(comp_l_trial),V_trial,output='tpd')
	# 			if tpd_v_trial<tpd_min_v:
	# 				tpd_min_v=tpd_v_trial
	# 				comp_v_min=comp_v_trial
	# 			if tpd_l_trial<tpd_min_l:
	# 				tpd_min_l=tpd_l_trial
	# 				comp_l_min=comp_l_trial
	# 			if delta_G_trial<delta_G_min:
	# 				delta_G_min=delta_G_trial
	# 			if tpd_min_l<0 or tpd_min_v<0 or delta_G_min<0:
	# 				second_counter=True
	# 				stable=False
	# 				c=0
	# 		else:
	# 			c2+=1
	# 		c+=1
	# 	if not stable:
	# 		V,comp_l,comp_h=V_trial,comp_v_trial,comp_l_trial
	# 		comps={'l':normalise_comp(comp_l),'h':normalise_comp(comp_h)}
			
	# 	else:
	# 		V,comps=V,{'h':comp}
	# 	stable=stable if not V in [0,1] else True

	# 	return (stable,V,comps)
	return stable,V,{'l':normalise_comp(comp_v),'h':normalise_comp(comp_l)}
def pvt_calculator(pvt,T,P):
	"""
	returns list of phase dicts from lightest to heaviest.
	phase dicts have the following keys:
	'type' (gas/liquid), 'comp' (dict of compositions),
	'props'(fluid properties),'fraction' (mole fraction),
	'quality' (mass fraction)
	{'comp':comps['l'],'props':{},'fraction':fraction,'quality':None,'type':'Unknown'}
	"""
	phases=[]
	#stable,fraction,comps=check_phase_stability(pvt,T,P)
	#stable,fraction,comps=michelson_stability(pvt,T,P)
	stable=False
	K_wilson=get_wilson_lnK(pvt,T,P)
	fraction,comp_v,comp_l=find_vapor_frcn(pvt,K_wilson)
	comps={'l':comp_v,'h':comp_l}
	if stable:
		if fraction<.5:
			phases.append({'comp':comps['h']})
			fraction=0.
		else:
			phases.append({'comp':comps['l']})
			fraction=1.
	else:
		N=len(pvt)
		tol=5e-4
		error=tol+1
		it=0
		comp_l,comp_h=comps['l'],comps['h']
		while error>tol and it<=200:
			try:
				##### newton #############
				comp_ln,comp_hn=flash_itterate(pvt,T,P,comp_l,comp_h,fraction)
			except:
				print "failed newton at P=%s, T=%s"%(P,T)
				try:
					fract,comp_ln,comp_hn=successive_substitution(pvt,T,P,comp_l,comp_h)
				except:
					print "failed successive substitution at P=%s, T=%s"%(P,T)
					break
			if it==20:
				print 'failed to converge at P=%s and T=%s' %(P,T)
			print comp_l,comp_h
			comp_ln,comp_hn=normalise_comp(comp_ln),normalise_comp(comp_hn)
			##### successive sub #####
			# try:
			# 	fract,comp_ln,comp_hn=successive_substitution(pvt,T,P,comp_l,comp_h)
			# 	comp_ln,comp_hn=normalise_comp(comp_ln),normalise_comp(comp_hn)
			# 	comp_l,comp_h=comp_ln,comp_hn
			# except:
			# 	print "failed successive substitution"
			# 	break
			# ##### new vapour fraction #####
			# lnK={}
			# for i in pvt:
			# 	lnK[i]=ln(comp_l[i])-ln(comp_h[i])
			# fraction=find_vapor_frcn(pvt,lnK)[0]
			error=max([abs(comp_ln[i]-comp_l[i]) for i in pvt]+[abs(comp_hn[i]-comp_h[i]) for i in pvt])
			comp_l,comp_h=normalise_comp(comp_ln),normalise_comp(comp_hn)
			it+=1
		phases.append({'comp':comp_l})
		phases.append({'comp':comp_h})
		K={}
		for i in pvt:
			K[i]=ln(comp_l[i]/comp_h[i])
		fraction=find_vapor_frcn(pvt,K)[0]
		print it
	##############
	return {'phases':phases,'V':fraction}
def michelson_stability(pvt,T,P,comp=None):
	if not comp:
		comp={}
		for i in pvt:
			comp[i]=pvt[i]['comp']
	lnfug_spec=fug_minimum_gibbs(pvt,T,P,1)
	lnK_wilson=get_wilson_lnK(pvt,T,P)
	### create vapour ###
	c,triv_v=0,False
	print "vapor:"
	while c<100 and not triv_v:
		lnK=lnK if not c==0 else lnK_wilson
		new_comp={}
		S_v=0
		for i in pvt:
			print '\t',comp[i],lnK[i]
			new_comp[i]=comp[i]*math.exp(lnK[i])
			S_v+=new_comp[i]
		new_comp=normalise_comp(new_comp)
		print '\t',new_comp
		lnfug_v=fug_minimum_gibbs(pvt,T,P,1,comp=new_comp)
		lnR={}
		for i in pvt:
			lnR[i]=ln(comp[i])+lnfug_spec[i]-ln(new_comp[i])-lnfug_v[i]-ln(S_v)
			lnK[i]=lnK[i]+lnR[i]
		if sum([lnK[i]**2 for i in pvt])<1e-4:
			triv_v=True
		if sum([(math.exp(lnR[i])-1)**2 for i in pvt])<1e-10:
			break
	comp_v=new_comp
	lnK_wilson=get_wilson_lnK(pvt,T,P)
	### create liquid ###
	c,triv_l=0,False
	print "liquid:"
	while c<100 and not triv_l:
		lnK=lnK if not c==0 else lnK_wilson
		new_comp_l={}
		S_l=0
		for i in pvt:
			print '\t',comp[i],lnK[i]
			new_comp_l[i]=comp[i]/math.exp(lnK[i])
			S_l+=new_comp_l[i]
		new_comp_l=normalise_comp(new_comp_l)
		print '\t',new_comp_l
		lnfug_l=fug_minimum_gibbs(pvt,T,P,1,comp=new_comp_l)
		lnR={}
		for i in pvt:
			lnR[i]=(ln(comp[i])+lnfug_spec[i]-ln(new_comp_l[i])-lnfug_l[i]-ln(S_l))
			lnK[i]=lnK[i]+lnR[i]
		if sum([lnK[i]**2 for i in pvt])<1e-4:
			triv_l=True
		if sum([(math.exp(lnR[i])-1)**2 for i in pvt])<1e-10:
			break
	comp_l=new_comp_l
	######## stability check #####
	if S_l<=1 and S_v<=1:
		stable=True
	elif triv_v and triv_l:
		stable=True
	elif triv_v and (S_l<=1):
		stable=True
	elif triv_l and (S_v<=1):
		stable=True
	else:
		stable=False
	return stable
def flash_itterate(pvt,T,P,comp_l,comp_h,fraction):
	""" 
	gauss-newton method for flash calc
	returns comp_l,comp_h
	"""
	comp={}
	for i in pvt:
		comp[i]=pvt[i]['comp']
	comp_l=normalise_comp(comp_l)
	comp_h=normalise_comp(comp_h)
	(lnfug_l,V_l),(lnfug_h,V_h)=fug_minimum_gibbs(pvt,T,P,fraction,comp=comp_l,returnV=True,phase='light'),fug_minimum_gibbs(pvt,T,P,(1- fraction),comp=comp_h,returnV=True,phase='heavy')
	################# Matrix solver ####################
	index={}		#maps index of variable to component name
	F,J,X,N={},{},{},len(pvt.keys())
	c=1
	for i in pvt.keys():
		index[c]=i
		c+=1
	c=1
	i_n=index[N-1]
	while c<=N:
		i=index[c]
		F_fug_dif=lnfug_l[i]+ln(comp_l[i])-lnfug_h[i]-ln(comp_h[i])
		F[c]=(F_fug_dif)
		if not c==N:
			X[c]=comp_l[i]
			X[c+N-1]=comp_h[i]
		if c<=N-2:
			F_mass_dif=(comp[i]- comp_l[i])/(comp_h[i]- comp_l[i])-(comp[i_n]- comp_l[i_n])/(comp_h[i_n]- comp_l[i_n])
			F[N+c]=F_mass_dif
		c+=1

	############ evaluate Jacobian ###########
	c,e=1,1e-8
	while c<=N:
		F_fug={}
		j=1
		name=index[c]
		while j<=N-1:
			# light
			comp_l1,comp_l2=update_dict(comp_l,index[j],comp_l[index[j]]-e),update_dict(comp_l,index[j],comp_l[index[j]]+e)
			lnfug_l1,lnfug_l2=fug_minimum_gibbs(pvt,T,P,fraction,comp=comp_l1,phase='light'),fug_minimum_gibbs(pvt,T,P,fraction,comp=comp_l2,phase='light')
			fug_difl1=lnfug_l1[name]+ln(comp_l1[name])-lnfug_h[name]-ln(comp_h[name])
			fug_difl2=lnfug_l2[name]+ln(comp_l2[name])-lnfug_h[name]-ln(comp_h[name])
			grad_l=(fug_difl2-fug_difl1)/(2*e)
			# heavy
			comp_h1,comp_h2=update_dict(comp_h,index[j],comp_h[index[j]]-e),update_dict(comp_h,index[j],comp_h[index[j]]+e)
			lnfug_h1,lnfug_h2=fug_minimum_gibbs(pvt,T,P,fraction,comp=comp_h1,phase='heavy'),fug_minimum_gibbs(pvt,T,P,fraction,comp=comp_h2,phase='heavy')
			fug_difh1=lnfug_l[name]+ln(comp_l[name])-lnfug_h1[name]-ln(comp_h1[name])
			fug_difh2=lnfug_l[name]+ln(comp_l[name])-lnfug_h2[name]-ln(comp_h2[name])
			grad_h=(fug_difh2-fug_difh1)/(2*e)
			F_fug[j]=grad_l
			F_fug[j+N-1]=grad_h
			#####################analytical###############
			# if j==c:
			# 	#light phase
			# 	try:
			# 		div_l_pos=ln_fug_div(pvt,T,P,V_l,comp_l,name)+1/comp_l[name]
			# 	except ZeroDivisionError:
			# 		div_l_pos=1e15
			# 	try:
			# 		div_l_neg=-ln_fug_div(pvt,T,P,V_h,comp_h,name)*fraction/(fraction-1)-fraction/(comp_h[name]*(1-fraction))
			# 	except ZeroDivisionError:
			# 		div_l_neg=-1e15
			# 	F_fug[j]=div_l_pos
			# 	#heavy phase
			# 	try:
			# 		div_h_neg=-ln_fug_div(pvt,T,P,V_h,comp_h,name)-1/comp_h[name]
			# 	except ZeroDivisionError:
			# 		div_h_neg=-1e15
			# 	try:
			# 		div_h_pos=ln_fug_div(pvt,T,P,V_l,comp_l,name)*(fraction-1)/fraction+(fraction-1)/(comp_l[name]*fraction)
			# 	except ZeroDivisionError:
			# 		div_h_pos=1e15
			# 	F_fug[j+N-1]=div_h_neg
			# else:
			# 	#light phase
			# 	try:
			# 		div_l_pos=ln_fug_div(pvt,T,P,V_h,comp_h,name)*fraction/(fraction-1)+fraction/(comp_h[name]*(1-fraction))
			# 	except ZeroDivisionError:
			# 		div_l_pos=1e15
			# 	try:
			# 		div_l_neg=-ln_fug_div(pvt,T,P,V_l,comp_l,name)-1/comp_l[name]
			# 	except ZeroDivisionError:
			# 		div_l_neg=-1e15
			# 	F_fug[j]=div_l_neg
			# 	#heavy phase
			# 	try:
			# 		div_h_pos=ln_fug_div(pvt,T,P,V_h,comp_h,name)+1/comp_h[name]
			# 	except ZeroDivisionError:
			# 		div_h_pos=1e15
			# 	try:
			# 		div_h_neg=-ln_fug_div(pvt,T,P,V_l,comp_l,name)*(fraction-1)/fraction-(fraction-1)/(comp_l[name]*fraction)
			# 	except ZeroDivisionError:
			# 		div_h_neg=-1e15
			# 	F_fug[j+N-1]=div_h_pos
			j+=1
		J[c]=F_fug
		if c<=N-2:
			F_mass={}
			j=1
			while j<=N-1:
				name_j=index[j]
				if j==c:
					#light
					div_l=(comp[name_j]-comp_h[name_j])/math.pow(comp_h[name_j]-comp_l[name_j],2)
					#heavy
					div_h=-(comp[name_j]-comp_l[name_j])/math.pow(comp_h[name_j]-comp_l[name_j],2)
					F_mass[j]=div_l
					F_mass[j+N-1]=div_h
				elif j==N-1:
					#light
					div_l=-(comp[name_j]-comp_h[name_j])/math.pow(comp_h[name_j]-comp_l[name_j],2)
					#heavy
					div_h=(comp[name_j]-comp_l[name_j])/math.pow(comp_h[name_j]-comp_l[name_j],2)
					F_mass[j]=div_l
					F_mass[j+N-1]=div_h
				else:
					F_mass[j]=0
					F_mass[j+N-1]=0
				j+=1
			J[c+N]=F_mass
		c+=1
	################## solve Newton-gauss #################
	print J,F
	print X
	J,F=gauss_elim(J,F)
	delta=solve_triangle_matrix(J,F)
	print delta
	comp_l={}
	comp_h={}
	i=1
	while i<=N-1:
		comp_l[index[i]]=-(delta[i]-X[i])
		if i==N-1:
			comp_l[index[N]]=1-sum([comp_l[j] for j in comp_l])
		i+=1
	while i<=2*N-2:
		comp_h[index[i-N+1]]=-(delta[i]-X[i])
		if i==2*N-2:
			comp_h[index[N]]=1-sum([comp_h[j] for j in comp_h])
		i+=1
	return comp_l,comp_h
def ln_fug_div(pvt,T,P,V,comp,i):
	eos_params=van_der_waal_coefs(pvt,T,P,comp=comp)
	v=V
	t=T+273.14
	d=8.314*t/P
	x=comp[i]
	a=eos_params['comp_params'][i]['a_i']
	m=eos_params['a_m']-a*x*x
	b=eos_params['comp_params'][i]['b_i']
	n=eos_params['b_m']-b*x
	div=-125*math.pow(3/2,2)*b*(b*(abs(a)*x*x+m)-2*(b*x+n)*math.sqrt(a*(abs(a)*x*x+m)))*ln((math.sqrt(2)*(b*x+n)+b*x+v+n)/(-math.sqrt(2)*(b*x+n)+b*x+v+n))/(4157*t*math.pow(b*x+n,3))+125*math.sqrt(2)*(-2*b*math.sqrt(a*(abs(a)*x*x+m))-2*a*abs(a)*x*(b*x+n)/math.sqrt(a*(abs(a)*x*x+m))+2*abs(a)*b*x)*ln((math.sqrt(2)*(b*x+n)+b*x+v+n)/(-math.sqrt(2)*(b*x+n)+b*x+v+n))/(4157*t*math.pow(b*x+n,2))-500*b*v*(abs(a)*x*x+m)/(4157*t*(2*math.pow(b*x+n,2)-math.pow(b*x+v+n,2)))+1000*abs(a)*v*x*(-b*x-n+b)/(4157*t*(2*math.pow(b*x+n,2)-math.pow(b*x+v+n,2)))-500*v*(-b*x-n+b)*(4*b*(b*x+n)-2*b*(b*x+v+n))*(abs(a)*x*x+m)/(4157*t*math.pow(2*math.pow(b*x+n,2)-math.pow(b*x+v+n,2),2))+125*math.sqrt(2)*(-math.sqrt(2)*(b*x+n)+b*x+v+n)*((math.sqrt(2)*b+b)/(-math.sqrt(2)*(b*x+n)+b*x+v+n)-(b-math.sqrt(2)*b)*(math.sqrt(2)*(b*x+n)+b*x+v+n)/math.pow(-math.sqrt(2)*(b*x+n)+b*x+v+n,2))*(b*(abs(a)*x*x+m)-2*(b*x+n)*math.sqrt(a*(abs(a)*x*x+m)))/(4157*t*math.pow(b*x+n,2)*(math.sqrt(2)*(b*x+n)+b*x+v+n))+b*(-b*x-n+b)/math.pow(-b*x+v-n,2)
	return div
def re_arrange_square_matrix(A,b):
	dimension=len(A)
	assert dimension==len(b), "matrix A and vector b have incompatible sizes"
	assert len(A[1])==dimension, "Not square matrix"
	graph={}
	row=1
	while row<=dimension:
	 	col=1
	 	U=[]
	 	while col<=dimension:
	 		if not A[row][col]==0:
	 			U.append(col)
	 		col+=1
	 	graph[row]=U
	 	row+=1
	matching,z,V=bipartite_match(graph)
	for i in range(1,dimension+1):
		assert i in V,"matrix for flash calx cant be solved. Check compositions"
	
	A_new,b_new={},{}
	for i in matching:
		A_new[i]=A[matching[i]]
		b_new[i]=b[matching[i]]
	return A_new,b_new
def bipartite_match(graph):
	'''Find maximum cardinality matching of a bipartite graph (U,V,E).
	The input format is a dictionary mapping members of U to a list
	of their neighbors in V.  The output is a triple (M,A,B) where M is a
	dictionary mapping members of V to their matches in U, A is the part
	of the maximum independent set in U, and B is the part of the MIS in V.
	The same object may occur in both U and V, and is treated as two
	distinct vertices if this happens.'''
	
	# initialize greedy matching (redundant, but faster than full search)
	matching = {}
	for u in graph: 
		for v in graph[u]:
			if v not in matching:
				matching[v] = u
				break
	
	while 1:
		# structure residual graph into layers
		# pred[u] gives the neighbor in the previous layer for u in U
		# preds[v] gives a list of neighbors in the previous layer for v in V
		# unmatched gives a list of unmatched vertices in final layer of V,
		# and is also used as a flag value for pred[u] when u is in the first layer
		preds = {}
		unmatched = []
		pred = dict([(u,unmatched) for u in graph])
		for v in matching:
			del pred[matching[v]]
		layer = list(pred)
		
		# repeatedly extend layering structure by another pair of layers
		while layer and not unmatched:
			newLayer = {}
			for u in layer:
				for v in graph[u]:
					if v not in preds:
						newLayer.setdefault(v,[]).append(u)
			layer = []
			for v in newLayer:
				preds[v] = newLayer[v]
				if v in matching:
					layer.append(matching[v])
					pred[matching[v]] = v
				else:
					unmatched.append(v)
		
		# did we finish layering without finding any alternating paths?
		if not unmatched:
			unlayered = {}
			for u in graph:
				for v in graph[u]:
					if v not in preds:
						unlayered[v] = None
			return (matching,list(pred),list(unlayered))

		# recursively search backward through layers to find alternating paths
		# recursion returns true if found path, false otherwise
		def recurse(v):
			if v in preds:
				L = preds[v]
				del preds[v]
				for u in L:
					if u in pred:
						pu = pred[u]
						del pred[u]
						if pu is unmatched or recurse(pu):
							matching[v] = u
							return 1
			return 0

		for v in unmatched: recurse(v)
def gauss_elim(A,b):
	dimension=len(A)
	row=1
	while row<=dimension:
		row2=row
		if A[row][row]==0:
			A_new,b_new={},{}
			row3=row
			while row3<=dimension:
				b_new[row3-row+1]=b[row3]
				col2=row
				row_new={}
				while col2<=dimension:
					row_new[col2-row+1]=A[row3][col2]
					col2+=1
				A_new[row3-row+1]=row_new
				row3+=1
			A_new,b_new=re_arrange_square_matrix(A_new,b_new)
			row3=row
			while row3<=dimension:
				col2=row
				b[row3]=b_new[row3-row+1]
				while col2<=dimension:
					A[row3][col2]=A_new[row3-row+1][col2-row+1]
					col2+=1
				row3+=1
		else:
			while row2<=dimension:
				row_dev={}
				if abs(A[row2][row])<1e-15:
					pass
				else:
					dev=float(A[row2][row])
					b[row2]=b[row2]/dev
					col=1
					while col<=dimension:
						row_dev[col]=A[row2][col]/dev
						col+=1
					if not row2==row:
						b[row2]=b[row2]-b[row]
						for i in row_dev:
							row_dev[i]=row_dev[i]-A[row][i]
					A[row2]=row_dev
				row2+=1
			row+=1
	return A,b
def solve_triangle_matrix(A,b):
	dimension=len(A)
	delta={}
	row=dimension
	while row>0:
		f=b[row]
		sum_f=0
		var=row
		while var<dimension:
			var+=1
			sum_f+=A[row][var]*delta[var]
		delta[row]=f-sum_f
		row-=1
	return delta
def update_dict(_dict,i,new_i):
	newDict={}
	for j in _dict:
		if i==j:
			newDict[j]=new_i
		else:
			newDict[j]=_dict[j]
	return newDict
def get_initial_bubble_point(pvt):
	P=1e5
	Ta=-273.1
	ga=rachford_rice_g(0,pvt,get_wilson_lnK(pvt,Ta,P))
	Tb=0
	gb=rachford_rice_g(0,pvt,get_wilson_lnK(pvt,Tb,P))
	while ga*gb>0 and not gb==0:
		Ta=Tb
		ga=gb
		Tb+=100
		gb=rachford_rice_g(0,pvt,get_wilson_lnK(pvt,Tb,P))
	if not gb==0:
		tol=1e-2
		e=tol+1
		while e>tol:
			Tc=(Ta+Tb)/2
			gc=rachford_rice_g(0,pvt,get_wilson_lnK(pvt,Tc,P))
			if gc*ga>0:
				Ta=Tc
				ga=gc
			else:
				Tb=Tc
				gb=gc
			e=abs(Ta-Tb)
		T=Tc
	else:
		T=Tb
	comp_v,comp_l={},{}
	K=get_wilson_lnK(pvt,T,P)
	return P,T,K
def itterate_envelope(pvt,P,T,lnK,spec,spec_value):
	disc,comp_l,comp_h=find_vapor_frcn(pvt,lnK)
	lnfug_l,lnfug_h=fug_minimum_gibbs(pvt,T,P,1,comp=comp_l,phase='light'),fug_minimum_gibbs(pvt,T,P,1,comp=comp_h,phase='heavy')
	################# Matrix solver ####################
	index={}		#maps index of variable to component name
	F,J,N={},{},len(pvt.keys())
	c=1
	for i in pvt.keys():
		index[c]=i
		if spec==i:
			spec=c
		c+=1
	########## fugacity continuity ############
	c=1
	e=1e-10
	while c<=N:
		i=index[c]
		F_fug_dif=lnfug_l[i]+ln(comp_l[i])-lnfug_h[i]-ln(comp_h[i])
		F[c]=F_fug_dif
		J[c]={}
		c2=1
		while c2<=N:
			j=index[c2]
			newK1,newK2=update_dict(lnK,j,lnK[j]-e),update_dict(lnK,j,lnK[j]+e)
			disc,comp_l1,comp_h1=find_vapor_frcn(pvt,newK1)
			disc,comp_l2,comp_h2=find_vapor_frcn(pvt,newK2)
			lnfug_l1,lnfug_h1=fug_minimum_gibbs(pvt,T,P,1,comp=comp_l1,phase='light'),fug_minimum_gibbs(pvt,T,P,1,comp=comp_h1,phase='heavy')
			lnfug_l2,lnfug_h2=fug_minimum_gibbs(pvt,T,P,1,comp=comp_l2,phase='light'),fug_minimum_gibbs(pvt,T,P,1,comp=comp_h2,phase='heavy')
			F_fug_dif1=lnfug_l1[i]+ln(comp_l1[i])-lnfug_h1[i]-ln(comp_h1[i])
			F_fug_dif2=lnfug_l2[i]+ln(comp_l2[i])-lnfug_h2[i]-ln(comp_h2[i])
			deriv_j=(F_fug_dif2- F_fug_dif1)/(2*e)
			J[c][c2]=deriv_j
			c2+=1
		########## vary P ##########
		P1,P2=math.exp(ln(P)-e),math.exp(ln(P)+e)
		lnfug_l1,lnfug_h1=fug_minimum_gibbs(pvt,T,P1,1,comp=comp_l,phase='light'),fug_minimum_gibbs(pvt,T,P1,1,comp=comp_h,phase='heavy')
		lnfug_l2,lnfug_h2=fug_minimum_gibbs(pvt,T,P2,1,comp=comp_l,phase='light'),fug_minimum_gibbs(pvt,T,P2,1,comp=comp_h,phase='heavy')
		F_fug_dif1=lnfug_l1[i]+ln(comp_l[i])-lnfug_h1[i]-ln(comp_h[i])
		F_fug_dif2=lnfug_l2[i]+ln(comp_l[i])-lnfug_h2[i]-ln(comp_h[i])
		deriv_p=(F_fug_dif2- F_fug_dif1)/(2*e)
		########## vary T ##########
		T1,T2=math.exp(ln(T+273.14)-e)-273.14,math.exp(ln(T+273.14)+e)-273.14
		lnfug_l1,lnfug_h1=fug_minimum_gibbs(pvt,T1,P,1,comp=comp_l,phase='light'),fug_minimum_gibbs(pvt,T1,P,1,comp=comp_h,phase='heavy')
		lnfug_l2,lnfug_h2=fug_minimum_gibbs(pvt,T2,P,1,comp=comp_l,phase='light'),fug_minimum_gibbs(pvt,T2,P,1,comp=comp_h,phase='heavy')
		F_fug_dif1=lnfug_l1[i]+ln(comp_l[i])-lnfug_h1[i]-ln(comp_h[i])
		F_fug_dif2=lnfug_l2[i]+ln(comp_l[i])-lnfug_h2[i]-ln(comp_h[i])
		deriv_t=(F_fug_dif2- F_fug_dif1)/(2*e)
		J[c][N+1]=deriv_p
		J[c][N+2]=deriv_t
		c+=1
	############### mass continuity ############
	F_mass_bal=sum([comp_l[i]-comp_h[i] for i in pvt])
	F[N+1]=F_mass_bal
	J[N+1]={}
	c2=1
	while c2<=N:
		j=index[c2]
		newK1,newK2=update_dict(lnK,j,lnK[j]-e),update_dict(lnK,j,lnK[j]+e)
		disc,comp_l1,comp_h1=find_vapor_frcn(pvt,newK1)
		disc,comp_l2,comp_h2=find_vapor_frcn(pvt,newK2)
		F_mass_bal1=sum([comp_l1[i]-comp_h1[i] for i in pvt])
		F_mass_bal2=sum([comp_l2[i]-comp_h2[i] for i in pvt])
		deriv_j=(F_mass_bal2- F_mass_bal1)/(2*e)
		J[N+1][c2]=deriv_j
		c2+=1
	J[N+1][N+1],J[N+1][N+2]=0,0
	######## specified variable ########
	J[N+2]={}
	if spec=="P":
		F[N+2]=ln(P)-ln(spec_value)
		c=1
		while c<=N+2:
			J[N+2][c]=0
			c+=1
		J[N+2][N+1]=-1
	elif spec=="T":
		F[N+2]=ln(T+273.14)-ln(spec_value+273.14)
		c=1
		while c<=N+2:
			J[N+2][c]=0
			c+=1
		J[N+2][N+2]=-1
	else:
		F[N+2]=lnK[index[spec]]-spec_value
		c=1
		while c<=N+2:
			J[N+2][c]=0
			c+=1
		J[N+2][spec]=-1
	################## solve Newton-gauss #################
	J,F=gauss_elim(J,F)
	delta=solve_triangle_matrix(J,F)
	lnK_new={}
	P_new=math.exp(ln(P)+delta[N+1])
	T_new=math.exp(ln(T+273.14)+delta[N+2])-273.14
	for i in index:
		lnK_new[index[i]]=delta[i]+lnK[index[i]]
	error=0
	delta_it={}
	for d in delta:
		error+=delta[d]**2
		if d<=N:
			delta_it[index[d]]=delta[d]
		elif d==N+1:
			delta_it['P']=delta[d]
		elif d==N+2:
			delta_it['T']=delta[d]
		else:
			assert False
	error=error**.5
	return P_new,T_new,lnK_new,error,delta_it
def calculate_envelope(pvt):
	P0,T0,lnK0=get_initial_bubble_point(pvt)
	tol=1e-7
	P,T,lnK,error,delta_it=itterate_envelope(pvt,P0,T0,lnK0,'P',P0)
	it=0
	delta_S={}
	while error>tol:
		P,T,lnK,error,delta_it=itterate_envelope(pvt,P,T,lnK,'P',P0)
		for i in pvt:
			delta_S[i]=delta_it[i]+ delta_S[i] if i in delta_S else delta_it[i]
		delta_S['P'],delta_S['T']=delta_it['P']+ delta_S['P'] if 'P' in delta_S else delta_it['P'],delta_it['T']+ delta_S['T'] if 'T' in delta_S else delta_it['T']
		it+=1
		print it,P,T,lnK,delta_S
	print delta_S
	return P,T,lnK,error,delta_it