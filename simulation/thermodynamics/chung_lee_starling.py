import math
exp=math.exp
sin=math.sin

k_i={'methanol':.215175,
   'ethanol':.174823,
   '1propanol':.143453,
   '2propanol':.143453,
   '1butanol':.131671,
   '2methyl1propanol':.131671,
   '1pentanol':.131671,
   '1hexanol':.121555,
   '1heptanol':.114230,
   'aceticacid':.091459,
   'ethanoicacid':.091459,
   'water':.075908}
mu_i={'methanol':1.7,
   'ethanol':1.69,
   '1propanol':1.68,
   '2propanol':1.66,
   '1butanol':1.66,
   '2methyl1propanol':1.64,
   '1pentanol':1.64,
   '1hexanol':1.63,
   '1heptanol':1.62,
   'aceticacid':1.74,
   'ethanoicacid':1.74,
   'water':1.85}
A,B,C,D,E,F,G,H,S,W=1.16145,.14874,.52487,.7732,2.16178,2.43787,-6.435e-4,7.27371,18.0323,-.7683
a0=[6.32402,.12102e-2,5.28346,6.62263,19.7454,-1.89992,24.2745,.79716,-.23816,.68629e-1]
a1=[50.4119,-.11536e-2,254.209,38.0957,7.63034,-12.5367,3.44945,1.11764,.67695e-1,.34793]
a2=[-51.6801,-.62571e-2,-168.481,-8.46414,-14.3544,4.98529,-11.2913,.12348e-1,-.8163,.59256]
a3=[1189.02,.37283e-1,3898.27,31.4178,31.5267,-18.1507,69.3466,-4.1161,4.02528,-.72663]
b0=[2.41657,-.50924,6.61069,14.5425,.79274,-5.8634,81.171]
b1=[.74824,-1.50936,5.62073,-8.91387,.82019,12.8005,114.158]
b2=[-.91858,-49.9912,64.7599,-5.63794,-.69369,9.58926,-60.841]
b3=[121.721,69.9834,27.0389,74.3435,6.31734,-65.5292,466.775]

def get_visc_and_conductivity(Phase):
    pvt,comp,T,mol_density=Phase.pvt,Phase.comp,Phase.T+273.15,Phase.mol_density
    rho=mol_density/1e6
    k=sum([sum([comp[i]*comp[j]*((k_i[i]*k_i[j])**.5 if (i in k_i and j in k_i) else 0) for j in comp]) for i in comp])
    sigma_i={i:.809*((pvt[i]['params']['v_c']*1000)**(1/3.)) for i in comp}
    sigma_ij={i:{j:(sigma_i[i]*sigma_i[j])**.5 for j in comp} for i in comp}
    sigma=(sum([sum([comp[i]*comp[j]*(sigma_ij[i][j]**3) for j in comp]) for i in comp]))**(1/3.)
    acc=sum([sum([comp[i]*comp[j]*.5*(pvt[i]['params']['acc']+pvt[j]['params']['acc'])*sigma_ij[i][j]**3 for j in comp]) for i in comp])/sigma**3
    epsilon_k_i={i:pvt[i]['params']['t_c']/1.2593 for i in comp}
    epsilon_k=sum([sum([comp[i]*comp[j]*((epsilon_k_i[i]*epsilon_k_i[j])**.5)*sigma_ij[i][j]**3 for j in comp]) for i in comp])/(sigma**3)
    M_ij={i:{j:1000*2*pvt[i]['params']['mm']*pvt[i]['params']['mm']/(pvt[i]['params']['mm']+pvt[i]['params']['mm']) for j in comp} for i in comp}
    M=((sum([sum([comp[i]*comp[j]*((epsilon_k_i[i]*epsilon_k_i[j])**.5)*(sigma_ij[i][j]**2)*(M_ij[i][j]**.5) for j in comp]) for i in comp]))/(epsilon_k*sigma**2))**2
    T_c=epsilon_k*1.2593
    V_c=(sigma/.809)**3
    T_s=T/epsilon_k
    mu=((sum([sum([comp[i]*comp[j]*((mu_i[i]**2)*(mu_i[j]**2) if (i in mu_i and j in mu_i) else 0)/((((epsilon_k_i[i]*epsilon_k_i[j])**.5)*(sigma_ij[i][j]**3))) for j in comp])for i in comp]))*(sigma**3)*epsilon_k)**.25
    mu_r=mu*131.3/((V_c*T_c)**.5)
    F_c=1-.2756*acc+.059035*(mu_r**4)+k
    Omega_s=A/(T_s**B)+C/exp(D*T_s)+E/exp(F*T_s)+G*(T_s**B)*sin(S*(T_s**W)-H)
    visc0=26.69e-6*F_c*((M*T)**.5)/((sigma**2)*(Omega_s))
    Y=rho*V_c/6.
    Ai=[a0[i]+a1[i]*acc+a2[i]*mu_r**4+a3[i]*k for i in range(10)]
    G1=(1.-.5*Y)/((1-Y)**3)
    G2=(Ai[0]*(1-exp(-Ai[3]*Y))/Y+Ai[1]*G1*exp(Ai[4]*Y)+Ai[2]*G1)/(Ai[0]*Ai[3]+Ai[1]+Ai[2])
    visc_p=(36.344e-6*((M*T_c)**.5)/(V_c**(2/3.)))*Ai[6]*(Y**2)*G2*exp(Ai[7]+Ai[8]/T_s+Ai[9]/(T_s**2))
    visc_k=visc0*(1/G2 +Ai[5]*Y)
    viscosity_mix=visc_k+visc_p
    ## conductivity ##
    alpha=(sum([comp[i]*pvt[i]['params']['heat_g']['C'](pvt[i]['params']['cg'], T-273.15)  for i in comp])/ 1e3 -8.314)/8.314
    beta=.7862-.7109*acc +1.3168*(acc**2)
    Z=2+10.5*((T/T_c)**2)
    Psi=1+alpha*((.215+.28288*alpha-1.061*beta+.26665*Z)/(.6366+beta*Z+1.061*alpha*beta))
    lambda0=7.452*(visc0/M)*Psi
    Bi = [b0[i] + b1[i] * acc + b2[i] * mu_r ** 4 + b3[i] * k for i in range(7)]
    H2=(Bi[0]*(1-exp(-Bi[3]*Y))/Y+Bi[1]*G1*exp(Bi[4]*Y)+Bi[2]*G1)/(Bi[0]*Bi[3]+Bi[1]+Bi[2])
    lambda_k=lambda0*((1/H2)+Bi[5]*Y)
    lambda_p=(3.039e-4*((T_c/M)**.5)/(V_c**(2/3.)))*Bi[6]*(Y**2)*H2*(T/T_c)**.5
    mixture_conductivity=lambda_k+lambda_p
    return (viscosity_mix*100,mixture_conductivity*4.1884*100)