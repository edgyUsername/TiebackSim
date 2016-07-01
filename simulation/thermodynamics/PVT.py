from simulation.thermodynamics import setup_comp_data as set_comp
import flash3
import pr_thermo

def get_props(pvt,T,P):
    sys=flash3.Flash(pvt,T,P)
    phases=sys.equilibriate()
    sys=pr_thermo.System(pvt,phases,T,P)
    props={'h':0}
    liq=[]
    c=0
    tot_mass=0
    m_l=0
    V_l=0
    Cp_l=0
    for p in sys.phases:
        if p.type=='gas':
            props['r_g']=p.mass_density
            props['m_g']=p.viscosity/1000
            props['Cp_g']=p.mass_Cp
        else:
            liq.append((p.mol_mass*p.mol_frac,p.viscosity/1000))
            m_l+=p.mol_mass*p.mol_frac
            V_l+=p.mol_frac/p.mol_density
            Cp_l+=p.mass_Cp*p.mol_mass*p.mol_frac
        tot_mass+=p.mol_mass*p.mol_frac
        props['h']+=p.mol_frac*p.mol_mass*p.mol_h
        c+=1
    props['Cp_l']=Cp_l/m_l
    props['r_l']=m_l/V_l
    props['vapQ']=1-m_l/tot_mass
    for i in range(len(liq)):
        if i==0:
            viscosity=liq[i][1]
            mass=liq[i][0]
        else:
            viscosity=viscosity*(2*viscosity+liq[i][1]-2*(viscosity-liq[i][1])*(liq[i][0]/sum([liq[i][0],mass])))/\
                      (2*viscosity+liq[i][1]+(viscosity-liq[i][1])*(liq[i][0]/sum([liq[i][0],mass])))
            mass+=liq[i][0]
    props['m_l']=viscosity

    return props