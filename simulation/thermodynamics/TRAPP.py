import math
R=8.314
ln=math.log
exp=math.exp
dens_r_c=.1628
pi=math.pi
sin=math.sin
def get_visc_and_conductivity(Phase):
    pvt,comp,density,mol_dens_i,phase=Phase.pvt,Phase.comp,Phase.mol_density,Phase.mol_dens_i,Phase.type
    T=Phase.T+273.15
    T_ri={i:(T)/pvt[i]['params']['t_c'] for i in comp}
    v_ri={i:1/(pvt[i]['params']['v_c']*1000*mol_dens_i[i]) for i in comp}
    w_i={i: pvt[i]['params']['acc'] for i in comp}
    z_ci={i: pvt[i]['params']['p_c']*1e5*pvt[i]['params']['v_c']/(1e3*R*pvt[i]['params']['t_c']) for i in comp}
    ############ Calculate psi ########################################################
    a2,b2,c2,d2=.394901,-1.023545,-.932813,-.754639
    v_ip={i:min(2,max(v_ri[i],.5)) for i in comp}
    T_ip={i:min(2,max(T_ri[i],.5)) for i in comp}
    G={i:a2*(v_ip[i]+b2)+c2*(v_ip[i]+d2)*ln(T_ip[i]) for i in comp}
    psi={i:(1+(w_i[i]-0.0115)*G[i])*0.2862130581857573/z_ci[i] for i in comp}
    ############ Calculate theta ######################################################
    a1,b1,c1,d1=.090569,-.862762,.316636,-.465684
    F={i:a1+b1*ln(T_ip[i])+(c1+d1/T_ip[i])*(v_ip[i]-.5) for i in comp}
    theta={i: 1+(w_i[i]-0.0115)*F[i] for i in comp}
    ############ Calculate h ##########################################################
    h_i={i:pvt[i]['params']['v_c']*psi[i]/0.0986 for i in comp}
    h_ij={i:{j:((h_i[j]**(1/3.)+h_i[i]**(1/3.))**3)/8. for j in comp} for i in comp}
    h=sum([sum([comp[i]*comp[j]*h_ij[i][j] for j in comp]) for i in comp])
    ############ Calculate f ##########################################################
    f_i={i: pvt[i]['params']['t_c']*theta[i]/190.564 for i in comp}
    f_ij={i:{j:(f_i[i]*f_i[j])**.5 for j in comp} for i in comp}
    f=sum([sum([comp[i]*comp[j]*f_ij[i][j]*h_ij[i][j] for j in comp]) for i in comp])/h
    ############ Calculate T_r and dens_r #############################################
    T_r=(T)/f
    dens_r=density*h*sum([comp[i]*pvt[i]['params']['mm']for i in comp])/1000
    ############ Calculate viscosity of reference fluid ###############################
    def _visc_ref(dens_r,T_r,M=.016):
        # dens_r=dens_r*M/1000
        a=[-10.238160427,174.22822961,17.460545674,-2847.6328289,.13368502192,142.07239767,5002.066972]
        b=[1.6969859271,-.13337234608,1.4,168.]
        c=[2.907741307e6,-3.312874033e6,1.608101838e6,-4.331904871e5,7.06248133e4,-7.11662075e3,4.3251744e2,-1.44591121e1,.2037119479]
        vr1=sum([c[i]*(T_r)**((i-3)/3.) for i in range(9)])
        vr2=b[0]+b[1]*(b[2]-ln((T_r)/b[3]))**2
        delta_vr=exp(a[0]+a[1]/(T_r))*(exp((a[2]+a[3]*(T_r)**(-1.5))*dens_r**.1
                                         +(dens_r/dens_r_c-1)*(dens_r**.5)*(a[4]+a[5]/(T_r)+a[6]*(T_r)**(-2.)))-1)
        visc_ref=vr1+vr2*dens_r+delta_vr
        return vr1,visc_ref #xE-4 cP
    vr1,visc_ref=_visc_ref(dens_r,T_r)
    ############ Calculate Translational conductivity for reference ###################
    a=[-7.197708227,85.67822264,12.471834689,-984.62522975,.35946850007,69.798412538,-872.88332851]
    b=[-.25276292,.33432859,1.12,168.]
    tcr1=15*R*vr1/(4*.016)
    tcr2=b[0]+b[1]*(b[2]-ln((T_r)/b[3]))**2
    delta_tcr=exp(a[0]+a[1]/(T_r))*(exp((a[2]+a[3]*(T_r)**(-1.5))*dens_r**.1
                                     +(dens_r/dens_r_c-1)*(dens_r**.5)*(a[4]+a[5]/(T_r)+a[6]*(T_r)**(-2.)))-1)
    trans_conductivity_ref=tcr1+tcr2*dens_r+delta_tcr
    ############ Calculate Viscosity for mixture ######################################
    M_ij = {i: {j: pvt[i]['params']['mm'] * pvt[j]['params']['mm'] * 2 / (pvt[i]['params']['mm'] + pvt[j]['params']['mm']) for j
        in comp} for i in comp}
    M_viscosity=sum([sum([comp[i]*comp[j]*(M_ij[i][j]**(.5))*(f_ij[i][j]**.5)*(h_ij[i][j]**(4/3.)) for j in comp]) for i in comp])**2/(f*h**(8/3.))
    F_viscosity=(M_viscosity*f/.016)**.5*h**(-2/3.)
    mixture_viscosity=visc_ref*F_viscosity

    ############ Calculate Translational conductivity for mixture #####################
    M_conductivity=f*h**(8/3.)*sum([sum([comp[i]*comp[j]*(M_ij[i][j]**(-.5))*(f_ij[i][j]**.5)*(h_ij[i][j]**(4/3.)) for j in comp]) for i in comp])**(-2)
    F_conductivity=(((.016/M_conductivity)*f)**.5)*h**(-2/3.)
    trans_conductivity=trans_conductivity_ref*F_conductivity
    ############ Calculate Internal conductivity for mixture ##########################
    Cp_ig={i: pvt[i]['params']['heat_g']['C'](pvt[i]['params']['cg'], (T-273.15))/1000 for i in comp}
    v_pure={i: _visc_ref(0,(T)/f_i[i])[1]*((pvt[i]['params']['mm']*f_i[i]/.016)**.5)*h_i[i]**(-4/3.-.5) for i in comp}
    c_i={i: 1.32*(Cp_ig[i]-5*R/2)*v_pure[i]/pvt[i]['params']['mm'] for i in comp}
    c_ij={i:{j: 2*c_i[i]*c_i[j]/(c_i[i]+c_i[j]) for j in comp} for i in comp}
    internal_conductivity=sum([sum([comp[i]*comp[j]*c_ij[i][j] for j in comp]) for i in comp])
    ############ Calculate Thermal conductivity for mixture ##########################
    mixture_conductivity=trans_conductivity+internal_conductivity # xE-7 W/m K
    ############ Calculate ENSKOG viscosity factor ##########################
    # M={i:{j:mol_dens_i[i]/(mol_dens_i[i]+mol_dens_i[j]) for j in comp} for i in comp}
    # sigma_i={i:(98.6*h_i[i]/(3.058e-10*6.023e23))**(1/3.) for i in comp}
    # sigma={i:{j:.5*(sigma_i[i]+sigma_i[j]) for j in comp}for i in comp}
    # b={i:{j:2*pi*(sigma[i][j]**3)/3. for j in comp} for i in comp}
    # g
    # y={i:{j:n*b[i][j]*g[i][j] for j in comp}for i in comp}
    # Y={i:comp[i]*(1+4*sum([comp[j]*M[j][i]*y[j][i] for j in comp])/5.) for i in comp}
    # ENSKOG_mix=sum([Beta[i]*Y[i] for i in ])+48*sum([sum([comp[i]*comp[j]*u[i][j]*y[i][j] for j in comp]) for i in comp])/(15*pi)
    return mixture_viscosity*1e-4,mixture_conductivity*1e-7
