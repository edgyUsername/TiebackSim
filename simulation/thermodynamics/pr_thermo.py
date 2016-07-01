import math
import peng_robinson
import TRAPP
import chung_lee_starling
R=8.314
ln= math.log
class Phase():
    def __init__(self,pvt,phase,T,P):
        self.comp=phase['comp']
        self.type=phase['type']
        self.pvt=pvt
        self.T=T
        self.P=P
        self.mol_frac=phase['frac']
        self._set_density()
        self._set_Cp()
        self._set_enthalpy()
        self._set_transport_props()

    def _set_density(self):
        pvt,T,P,comp=self.pvt,self.T,self.P,self.comp
        # gas
        Z, self.eos_params = peng_robinson.solve_PR_for_Z(pvt, T, P, comp=comp)
        Z = Z[-1]
        self.g_mol_density = P / (Z * R * (T + 273.15))
        self.mol_mass = sum([comp[i] * pvt[i]['params']['mm'] for i in comp])
        # liquid
        self.l_mol_density = 1e3 / sum \
            ([comp[i] * (pvt[i]['params']['v_c'] * (.29056 - .08775 * pvt[i]['params']['acc']) ** (
            (1 - (T + 273.15) / pvt[i]['params']['t_c']) * abs(1 - (T + 273.15) / pvt[i]['params']['t_c']) ** (
            2 / 7. - 1))) for i in comp])
        if self.type=='gas':
            self.mol_density=self.g_mol_density
            self.mol_dens_i={i:P/(peng_robinson.solve_PR_for_Z({i:pvt[i]},T,P,comp={i:1.})[0][-1]*R*(T+273.15)) for i in comp}
            self.mass_density=self.mol_mass*self.g_mol_density
        elif self.type in ['hc','aq']:
            self.mol_density=self.l_mol_density
            self.mol_dens_i={i:1e3/comp[i]*(pvt[i]['params']['v_c']*(.29056-.08775*pvt[i]['params']['acc'])**((1-(T+273.15)/pvt[i]['params']['t_c'])*abs(1-(T+273.15)/pvt[i]['params']['t_c'])**(2/7.-1))) for i in comp}
            self.mass_density=self.mol_mass*self.l_mol_density
        else:
            assert False
    def _set_Cp(self):
        pvt, comp, T = self.pvt, self.comp, self.T
        self.mol_Cp=sum([comp[i]*pvt[i]['params']['heat_l']['C'](pvt[i]['params']['cl'],T,pvt[i]['params']['t_c']) for i in comp])/1e3\
                if self.type in ['hc','aq'] else sum([comp[i]*pvt[i]['params']['heat_g']['C'](pvt[i]['params']['cg'], T)  for i in comp])/ 1e3
        self.mass_Cp=self.mol_Cp/self.mol_mass

    def _set_enthalpy(self):
        pvt,comp,T,P=self.pvt,self.comp,self.T,self.P
        h_i_g=sum([1e-3*comp[i] * pvt[i]['params']['heat_g']['h'](pvt[i]['params']['cg'], T) for i in comp])\
              +sum([comp[i]*pvt[i]['params']['hf']])
        v=max(1/self.mol_density,self.eos_params['b_m']+1e-16)
        vg=1/self.g_mol_density
        vl=max(1/self.l_mol_density,self.eos_params['b_m']+1e-16)
        # v_ig=R*(T+273.15)/P
        a,b=[self.eos_params[i] for i in ['a_m','b_m']]
        h_dep=.5*a*ln(v**2+2*b*v-b**2)/(v+b)+R*(T+273.15)*(2*ln(v/(v-b))+(P*v/(R*(T+273.15)))-1)
        h=h_i_g+h_dep
        self.u_vap=(.5*a*ln(vg**2+2*b*vg-b**2)/(vg+b)+R*(T+273.15)*(2*ln(vg/(vg-b))+(P*vg/(R*(T+273.15)))-1))\
                   -(.5*a*ln(vl**2+2*b*vl-b**2)/(vl+b)+R*(T+273.15)*(2*ln(vl/(vl-b))+(P*vl/(R*(T+273.15)))-1))\
                   -P*(vg-vl)
        self.mol_h=h
        self.mass_h=h/self.mol_mass
    def _set_transport_props(self):
        # visc in cP and cond in W/(m K)
        if self.type=='gas':
            self.viscosity, self.thermal_conductivity =TRAPP.get_visc_and_conductivity(self)
        else:
            self.viscosity,self.thermal_conductivity=chung_lee_starling.get_visc_and_conductivity(self)

class System():
    def __init__(self,pvt,phases,T,P):
        self.pvt=pvt
        self.T=T
        self.P=P
        self.phases=[Phase(pvt,p,T,P) for p in phases]