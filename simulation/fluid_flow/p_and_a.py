import sys
import math
from os import path
sys.path.append( path.dirname( path.dirname(path.dirname( path.abspath(__file__)) ) ) )
#import settings
from simulation.fluid_flow import g_s

def estimateSigma (r_l,T,P):
	T_F=T*1.8+32
	if T_F<68:
		T_F=68
	elif T_F>100:
		T_F=100
	sg=r_l/999.016
	API=141.5/sg-131.5
	sigma68=39-0.2571*API
	sigma100=37.5-0.2571*API
	sigma=sigma68-(T_F-68)*(sigma68-sigma100)/32
	p=P/6894.757
	sigma=sigma*math.exp(-8.6306e-4*p)
	return sigma/1000

def itterate_h(h,r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta):
	### derivation ###
	# -Al(dp/dL)-twlSl+tiSi-r_l Al g/gc sin theta=
	# -Ag(dp/dL)-twgSg-tiSi-r_g Ag g/gc sin theta
	# -(dp/dL)-twl Sl/Al+ti Si/Al-r_l g/gc sin theta=0
	# -twl*Sl/(Al*twg)+ Sg/Ag+ ti*Si*(1/Al+1/Ag)/twg+(r_g-r_l)*g*sin(theta)/twg=0
	d=ID
	Al=d*d*0.25*(math.pi- math.acos(2*h-1)+(2*h-1)*math.sqrt(1-(2*h-1)*(2*h-1)))
	Ag=d*d*0.25*(math.acos(2*h-1)-(2*h-1)*math.sqrt(1-(2*h-1)*(2*h-1)))
	Sl=d*(math.pi- math.acos(2*h-1))
	Sg=d*math.acos(2*h-1)
	Si=d*math.sqrt(1-(2*h-1)*(2*h-1))
	dg=4*Ag/(Sg+Si)
	H_is=Al/(Al+Ag)

	Re_g=dg*r_g*(v_sg/(1-H_is))/m_g
	f_g=g_s.calcDWf(roughness,ID,Re_g)
	t_wG=f_g*r_g*math.pow(v_sg/(1-H_is),2)/2

	Re_sl=ID*r_l*v_sl/m_l
	f_sl=g_s.calcDWf(roughness,ID,Re_sl)
	f_l=.452*math.pow(f_sl,0.731)
	t_wL=f_l*r_l*math.pow(v_sl/H_is,2)/2

	Vi=v_sg/(1-H_is)-v_sl/H_is
	Fr=v_sl/(H_is*math.sqrt(9.81*h*d))
	f_i=(0.004+0.5e-6*Re_sl)*math.pow(Fr,1.335)*r_l*9.81*ID/(r_g*math.pow(v_sg/(1-H_is),2))
	t_i=f_i*r_g*Vi*abs(Vi)/2

	y=-t_wL*Sl/(Al*t_wG)+ Sg/Ag+ t_i*Si*(1/Al+1/Ag)/t_wG+(r_g-r_l)*9.81*math.sin(theta)/t_wG
	return (y,Al,Ag,Sl,Sg,Si,t_wL,t_wG,t_i)

def calc_h(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta):
	### illinois regula falsi - iterate h###
	## edited to bisection ##
	tol=0.001
	e=tol+1
	params=(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta)
	a=0.00001
	it_a=itterate_h(a,*params)
	f_a=it_a[0]

	b=.99999

	b=a+0.001
	it_b=itterate_h(b,*params)
	f_b=it_b[0]
	while f_b*f_a>=0 and b<.99999:
		b+=.001
		it_b=itterate_h(b,*params)
		f_b=it_b[0]

	while e>tol:
		c=(a+b)/2
		it_c=itterate_h(c,*params)
		f_c=it_c[0]

		if (f_a<0) == (f_c<0):
			a=c
			f_a=f_c
		else:
			b=c
			f_b=f_c
		e=abs((a-b)/c)
	y,Al,Ag,Sl,Sg,Si,t_wL,t_wG,t_i=it_c
	h=c
	return (h,Al,Ag,Sl,Sg,Si,t_wL,t_wG,t_i)

def itterate_d(d,r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P):
	area=math.pi*math.pow(ID/2,2)
	Af=math.pi*d*ID*(ID-d*ID)
	Ac=math.pi*math.pow(ID-2*d*ID,2)/4
	Sl=math.pi*ID
	Si=math.pi*ID*(1-2*d)

	sigma=estimateSigma(r_l,T,P)
	NB=math.pow(m_l*v_sg/sigma,2)*r_g/r_l
	FE=1/(1+1/(0.735*math.pow(NB,0.074)*math.pow(v_sg/v_sl,0.2)))
	
	D_f=4*d*(ID-d*ID)
	v_f=v_sl*area*(1-FE)/Af
	Re_f=D_f*r_l*v_f/m_l
	f_f=g_s.calcDWf(roughness,ID,Re_f)
	t_wL=f_f*r_l*math.pow(v_f,2)/2
	
	v_c=(v_sg+v_sl*FE)*area/Ac
	mass_l=v_sl*area*r_l
	mass_v=v_sg*area*r_g
	r_c=(FE*mass_l+mass_v)/(FE*mass_l/r_l+mass_v/r_g)
	Re_c=(ID-2*d*ID)*v_c*r_c/m_g
	f_c=g_s.calcDWf(0,ID-2*d*ID,Re_c)
	f_i=0.24*math.pow(sigma/(r_c*(ID-2*d*ID)*math.pow(v_c,2)),0.085)*math.pow(Re_f,0.305)

	t_i=f_i*r_c*(v_c-v_f)*abs(v_c-v_f)/2
	y=t_wL*Sl/Af-t_i*Si*(1/Af+1/Ac)+(r_l-r_c)*9.81*math.sin(theta)

	return (y,Af,Ac,Sl,Si,t_wL,t_i,FE)

def calc_d(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P):
	### illinois regula falsi - iterate h###
	## edited to bisection ##
	tol=0.001
	e=tol+1
	params=(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P)
	a=0.00001
	it_a=itterate_d(a,*params)
	f_a=it_a[0]
	
	b=a+0.001
	it_b=itterate_d(b,*params)
	f_b=it_b[0]

	while f_b*f_a>=0 and b<.495:
		b+=.001
		it_b=itterate_d(b,*params)
		f_b=it_b[0]

	assert f_a*f_b<0, "film thickness convergance error"
	while e>tol:
		c=(a+b)/2
		it_c=itterate_d(c,*params)
		f_c=it_c[0]

		if (f_a<0) == (f_c<0):
			a=c
			f_a=f_c
		else:
			b=c
			f_b=f_c
		e=abs((a-b)/c)
	y,Af,Ac,Sl,Si,t_wL,t_i,FE=it_c
	d=c
	return (d,Af,Ac,Sl,Si,t_wL,t_i,FE)

def itterate_d_min(d,r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P):
	area=math.pi*math.pow(ID/2,2)
	Af=math.pi*d*ID*(ID-d*ID)
	Ac=math.pi*math.pow(ID-2*d*ID,2)/4
	Sl=math.pi*ID
	Si=math.pi*ID*(1-2*d)
	Vf=Af/area

	sigma=estimateSigma(r_l,T,P)
	NB=math.pow(m_l*v_sg/sigma,2)*r_g/r_l
	FE=1/(1+1/(0.735*math.pow(NB,0.074)*math.pow(v_sg/v_sl,0.2)))
	v_c=(v_sg+v_sl*FE)*area/Ac
	mass_l=v_sl*area*r_l
	mass_v=v_sg*area*r_g
	r_c=(FE*mass_l+mass_v)/(FE*mass_l/r_l+mass_v/r_g)

	D_f=4*d*(ID-d*ID)
	v_f=v_sl*area*(1-FE)/Af
	Re_f=D_f*r_l*v_f/m_l
	f_f=g_s.calcDWf(roughness,ID,Re_f)
	
	y=2*f_f*r_l*math.pow(v_sl*(1-FE),2)/((r_l-r_c)*9.81*ID*math.sin(theta))-math.pow(Vf,3)*(1-1.5*Vf)/(2-1.5*Vf)

	return y

def calc_d_min(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P):
	### illinois regula falsi - iterate h###
	## edited to bisection ##
	tol=0.001
	e=tol+1
	params=(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P)
	a=0.000001
	f_a=itterate_d_min(a,*params)

	b=a+0.001
	f_b=itterate_d_min(b,*params)
	while f_b*f_a>=0 and b<.497:
		b+=.001
		f_b=itterate_d_min(b,*params)

	print f_a,f_b
	assert f_a*f_b<0, "min film thickness convergance error"
	while e>tol:
		c=(a+b)/2
		f_c=itterate_d_min(c,*params)

		if (f_a<0) == (f_c<0):
			a=c
			f_a=f_c
		else:
			b=c
			f_b=f_c
		e=abs((a-b)/c)
	d=c
	return d

def calcPressureDrop(P0,m,vapQ,r_g,r_l,m_g,m_l,z,roughness,ID,L,T):
	bp={"C_lift":.8,"gamma":1.3,"db":7e-3} #C_lift=.4 to 1.2;gamma=1.1 to 1.5, db= 4-10 mm
	par=[float(i) for i in [P0,m,vapQ,r_g,r_l,m_g,m_l,z,roughness,ID,L,T]]
	P0,m,vapQ,r_g,r_l,m_g,m_l,z,roughness,ID,L,T=par

	Q_l=m*(1-vapQ)/r_l
	Q_g=m*vapQ/r_g
	H=Q_l/(Q_g+Q_l)
	A=math.pi*math.pow(ID,2)/4
	v_sg=Q_g/A
	v_sl=Q_l/A
	v_m=v_sg+v_sl
	theta=math.asin(z/L)
	regime=None

	sigma=estimateSigma(r_l,T,P0)
	# check "dispersed bubble"
	E_ls=1/(1+math.pow(v_m/8.66,1.39))
	C_g=v_sg/v_m
	if E_ls<0.48 and C_g<=0.52:
		regime="dispersed bubble"
	else:
		# check "stratified"
		if z<=0:
			h,Al,Ag,Sl,Sg,Si,t_wL,t_wG,t_i=calc_h(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta)
			cos_t=math.cos(theta)
			cos_t=cos_t if cos_t>0.02 else 0.02
			dAl=-4*h*ID*(h*ID-ID)/(ID*math.sqrt(1-math.pow(2*h-1,2)))
			v_bridge=(1-h*ID)*math.sqrt((r_l-r_g)*9.81*Ag*cos_t/(r_g*dAl))
			H_is=Al/(Al+Ag)
			v_g=v_sg/(1-H_is)
			v_l=v_sl/H_is

			Re_sl=ID*r_l*v_sl/m_l
			f_sl=g_s.calcDWf(roughness,ID,Re_sl)
			f_l=.452*math.pow(f_sl,0.731)

			if v_g<=v_bridge and v_l<=math.sqrt(9.81*ID*(1-h)*cos_t/f_l):
				if v_g<=math.sqrt(4*m_l*(r_l-r_g)*9.81*cos_t/(0.06*r_l*r_g*v_l)) and v_l/math.sqrt(9.81*h*ID)<=1.4:
					regime="stratified smooth"
				else:
					regime="stratified wavy"
		if regime==None:
			# check "annular mist"
			try:
				d,Af,Ac,Sl,Si,t_wL,t_i,FE=calc_d(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P0)
				d_min=calc_d_min(r_l,r_g,m_l,m_g,v_sl,v_sg,roughness,ID,theta,T,P0)
				v_c=(v_sg+v_sl*FE)*(Ac+Af)/Ac
				v_f=v_sl*(Ac+Af)*(1-FE)/Af
				E_l_mist=v_sl*(v_f*FE+v_c*(1-FE))/(v_c*v_f)
				if d<d_min and E_l_mist<=.24:
					regime="annular mist"
			except AssertionError:
				print "film calc conversion error- skipping"
			if regime==None and E_ls<.48:			
				# check "froth"
				regime="froth"
			elif regime==None:
				VB=1.53*math.pow(9.81*(r_l-r_g)*sigma/(r_l*r_l),.25)*math.sin(theta)
				Re_mL=r_l*v_m*ID/m_l
				C0=(1.64+.12*math.sin(theta))*math.pow(Re_mL,-.031)
				VGdb=C0*v_m+VB

				Bo=(r_l-r_g)*9.81*ID*ID/sigma
				beta=Bo*math.exp(3.278-1.424*math.log(Bo,math.e))
				Vdv_inf=0.345*(1-math.exp(-beta))*math.sqrt(9.81*ID*(r_l-r_g)/r_l)
				Vdh_inf=(.54-1.76/math.pow(Bo,.56))*math.sqrt(9.81*ID*(r_l-r_g)/r_l)
				Vd_inf=Vdh_inf*math.cos(theta)+Vdv_inf*math.sin(theta)
				Re_inf=r_l*Vd_inf*ID/(2*m_l)
				f_m=.316*math.sqrt(Re_inf)
				f_m=f_m if f_m<1 else 1
				Vd=Vd_inf*f_m
				Vt=C0*v_m+Vd
				E_l_slug=(E_ls*Vt+VGdb*(1-E_ls)-v_sg)/Vt
				E_l_slug=E_l_slug if E_l_slug<=1 else 1-C_g
				# check "bubble"
				v_bubble=1.41*math.pow(9.81*(r_l-r_g)*sigma/(r_l*r_l),.25)*math.sin(theta)
				cond_2=math.cos(theta)<=(3/(4*math.sqrt(2)))*v_bubble*v_bubble*bp['C_lift']*bp['gamma']*bp['gamma']/(9.81*bp['db'])
				cond_1=ID>19*math.sqrt((r_l-r_g)*sigma/(r_l*r_l*9.81))
				cond_3=E_l_slug>.75
				if cond_3 and cond_1 and cond_2:
					regime="bubble"
				elif E_l_slug<=.24:
					regime="froth"
				elif E_ls<=.9:
					regime="slug"
				else:
					regime="elongated bubble"
	print regime


