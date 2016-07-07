import math
import peng_robinson
import numpy as np
from scipy.optimize import minimize_scalar,minimize
aq_comps=['water','ethanol','methanol','carbondioxide']
R=8.314
ln=math.log
exp=math.exp
def generate_pvt(pvt,components):
    """
    given a composition, returns the new pvt object with updated normalized composition
    :param pvt: pvt dict
    :param components: list of components
    :return: pvt dict
    """
    new_pvt={}
    total_moles=0
    for i in components:
        total_moles+=pvt[i]['comp']
    for i in components:
        pvt_i=pvt[i].copy()
        pvt_i['comp']=pvt_i['comp']/total_moles
        new_pvt[i]=pvt_i
    return new_pvt
def normalize_comp(comp):
    tot=sum([comp[i] for i in comp])
    return {i: comp[i]/tot for i in comp}


class Flash():

    def __init__(self,pvt,T,P):
        self.pvt=pvt
        self.T=T
        self.P=P
        self.phases=[]
        self.index=[i for i in pvt]

    def _split(self):
        lnK=peng_robinson.get_wilson_lnK(self.pvt,self.T,self.P)
        V,comp_v,comp_l=peng_robinson.find_vapor_frcn(self.pvt,lnK)
        self.phases=[{'frac':V,'comp':comp_v},{'frac':1-V,'comp':comp_l}]
    def _successive_sub(self):
        fug_v=peng_robinson.fug_minimum_gibbs(self.pvt,self.T,self.P,1.,comp=normalize_comp(self.phases[0]['comp']),phase='light')
        fug_l=peng_robinson.fug_minimum_gibbs(self.pvt,self.T,self.P,1.,comp=normalize_comp(self.phases[1]['comp']),phase='heavy')
        lnK={i:fug_l[i]-fug_v[i] for i in self.pvt}
        V,comp_v,comp_l=peng_robinson.find_vapor_frcn(self.pvt,lnK)
        self.phases=[{'frac':V,'comp':comp_v},{'frac':1-V,'comp':comp_l}]
        return ([(exp(fug_l[i])*comp_l[i]/(exp(fug_v[i])*comp_v[i])) for i in self.index],lnK)

    def _a_successive_sub(self,R_old,K_old):
        fug_v = peng_robinson.fug_minimum_gibbs(self.pvt, self.T, self.P, 1., comp=self.phases[0]['comp'],
                                                phase='light')
        fug_l = peng_robinson.fug_minimum_gibbs(self.pvt, self.T, self.P, 1., comp=self.phases[1]['comp'],
                                                phase='heavy')
        lnK_new = {i: fug_l[i] - fug_v[i] for i in self.pvt}
        V, comp_v, comp_l = peng_robinson.find_vapor_frcn(self.pvt, lnK_new)
        R_new=[(exp(fug_l[i])*comp_l[i]/(exp(fug_v[i])*comp_v[i])) for i in self.index]
        lambda_i=[(R_old[i]-1)/(R_old[i]-R_new[i]) for i in range(len(self.index))]
        lnK_new={self.index[i]:K_old[self.index[i]]*R_old[i]**lambda_i[i] for i in range(len(self.index))}
        V, comp_v, comp_l = peng_robinson.find_vapor_frcn(self.pvt, lnK_new)
        self.phases = [{'frac': V, 'comp': comp_v}, {'frac': 1 - V, 'comp': comp_l}]
        return R_new,lnK_new

    def _tangent_plane_analysis(self):
        lnK=peng_robinson.get_wilson_lnK(self.pvt,self.T,self.P)
        comp_l={i:self.pvt[i]['comp']/exp(lnK[i]) for i in self.pvt}
        comp_v={i:self.pvt[i]['comp']*exp(lnK[i]) for i in self.pvt}
        comp_l=normalize_comp(comp_l)
        comp_v=normalize_comp(comp_v)
        fug_spec_l=peng_robinson.fug_minimum_gibbs(self.pvt,self.T,self.P,1.,phase='light')
        fug_spec_h=peng_robinson.fug_minimum_gibbs(self.pvt,self.T,self.P,1.,phase='heavy')
        z={i:self.pvt[i]['comp'] for i in self.pvt}
        fug_spec=fug_spec_h #if sum([z[i]*(fug_spec_h[i]-fug_spec_l[i]) for i in z]) <0 else fug_spec_l
        tpd_min_l=None
        x_min_l=None
        triv_l=False
        for it in range(100):
            fug_l=peng_robinson.fug_minimum_gibbs(self.pvt,self.T,self.P,1.,comp=comp_l,phase='heavy')
            tpd_l=sum([comp_l[i]*(ln(comp_l[i])+fug_l[i]-fug_spec[i]-ln(z[i])) for i in self.pvt])
            if tpd_min_l ==None:
                tpd_min_l=tpd_l
                x_min_l=comp_l
            else:
                if tpd_l<tpd_min_l:
                    tpd_min_l = tpd_l
                    x_min_l = comp_l
            comp_l_n={i:z[i]*exp(fug_spec[i]-fug_l[i]) for i in z}
            comp_l_n=normalize_comp(comp_l_n)
            if tpd_min_l < 0:
                for it2 in range(3):
                    fug_l = peng_robinson.fug_minimum_gibbs(self.pvt, self.T, self.P, 1., comp=comp_l, phase='heavy')
                    comp_l = {i: z[i] * exp(fug_spec[i] - fug_l[i]) for i in z}
                    comp_l=normalize_comp(comp_l)
                break
            elif sum([(comp_l[i] - comp_l_n[i]) ** 2 for i in z]) < 1e-6:
                break
            elif sum([(comp_l_n[i] - z[i]) ** 2 for i in z]) < 1e-6:
                triv_l = True
                break
            comp_l = comp_l_n.copy()
        tpd_min_v = None
        x_min_v = None
        triv_v=False
        for it in range(100):
            fug_v = peng_robinson.fug_minimum_gibbs(self.pvt, self.T, self.P, 1., comp=comp_v, phase='light')
            tpd_v=sum([comp_v[i]*(ln(comp_v[i])+fug_v[i]-fug_spec[i]-ln(z[i])) for i in self.pvt])
            if tpd_min_v ==None:
                tpd_min_v=tpd_v
                x_min_v=comp_v
            else:
                if tpd_v<tpd_min_v:
                    tpd_min_v = tpd_v
                    x_min_v = comp_v
            comp_v_n={i:z[i]*exp(fug_spec[i]-fug_v[i]) for i in z}
            comp_v_n=normalize_comp(comp_v_n)
            if tpd_min_v<0:
                for it2 in range(3):
                    fug_v = peng_robinson.fug_minimum_gibbs(self.pvt, self.T, self.P, 1., comp=comp_v, phase='light')
                    comp_v={i:z[i]*exp(fug_spec[i]-fug_v[i]) for i in z}
                    comp_v=normalize_comp(comp_v)
                break
            elif sum([(comp_v[i]-comp_v_n[i])**2 for i in z])<1e-6:
                break
            elif sum([(comp_v_n[i]-z[i])**2 for i in z])<1e-6:
                triv_v=True
                break
            comp_v=comp_v_n.copy()

    def flash(self):
        index=[i for i in self.pvt]
        z=[self.pvt[i]['comp'] for i in index]
        n=len(index)
        x0=[ln(self.phases[0]['comp'][i]) for i in index]+[ln(self.phases[1]['comp'][i]) for i in index]
        def fun(x):
            comp_l={index[i]:exp(x[n+i]) for i in range(n)}
            comp_v={index[i]:exp(x[i]) for i in range(n)}
            fug_v=peng_robinson.fug_minimum_gibbs(self.pvt,self.T,self.P,1.,comp=comp_v,phase='light')
            fug_l=peng_robinson.fug_minimum_gibbs(self.pvt,self.T,self.P,1.,comp=comp_l,phase='heavy')
            F=[fug_v[i]+ln(comp_v[i])-fug_l[i]-ln(comp_l[i]) for i in index]
            F+=[(z[i]-exp(x[i]))/(exp(x[n+i])-exp(x[i]))-(z[n-1]-exp(x[n-1]))/(exp(x[-1])-exp(x[n-1])) for i in range(n-1)]
            return np.dot(F,F)

        cons=({'type': 'ineq',
               'fun': lambda x: np.array([-c for c in x]),
               'jac': lambda x: np.array([[-1. if c1 == c2 else 0 for c1 in range(len(x))] for c2 in range(len(x))])},
              {'type': 'eq',
               'fun': lambda x: np.array([1-sum([exp(c) for c in x[:n]]),1-sum([exp(c) for c in x[n:]])]),
               'jac': lambda x: np.array([[-exp(c) if c in x[:n] else 0 for c in x],[-exp(c) if c in x[n:] else 0 for c in x]])
               })
        res=minimize(fun,x0,method='SLSQP',constraints=cons)
        comp_v = [exp(c) for c in res.x[:n]]
        comp_l = [exp(c) for c in res.x[n:]]

        def func(frac):
            return ((1 - frac) ** 2) * np.dot(comp_l, comp_l) + 2 * frac * (1 - frac) * np.dot(comp_l, comp_v) +\
                   (frac ** 2) * np.dot(comp_v, comp_v) - np.dot(z, z) \
                   - 2 * (1 - frac) * np.dot(z, comp_l) - 2 * frac * np.dot(z, comp_v)
        bounds = [0., 1.]
        fraction = minimize_scalar(func, bounds=bounds, method='Bounded').x
        comp_v={index[i]:comp_v[i] for i in range(n)}
        comp_l={index[i]:comp_l[i] for i in range(n)}
        self.phases = [{'frac': fraction, 'comp': comp_v}, {'frac': 1 - fraction, 'comp': comp_l}]
    def _split_aq(self):
        gas=self.phases[0]
        liq=self.phases[1]
        aq_comp={}
        hc_comp={}
        aq_frac=0
        for i in liq['comp']:
            if i in aq_comps:
                aq_comp[i]=liq['comp'][i]
                aq_frac+=liq['comp'][i]
            else:
                hc_comp[i]=liq['comp'][i]
        self.phases=[{'type':'gas','frac':gas['frac'],'comp':gas['comp']},
                     {'type': 'hc', 'frac': (1-aq_frac)*liq['frac'], 'comp': normalize_comp(hc_comp)},
                     {'type': 'aq', 'frac': aq_frac* liq['frac'], 'comp': normalize_comp(aq_comp)}]
        for p in range(len(self.phases)):
            if self.phases[p]['frac']==0:
                self.phases.pop(p)

    def equilibriate(self):
        self._split()
        R_old,K_old=self._successive_sub()
        frac_old = self.phases[0]['frac']
        tol=1e-8
        error=tol+1
        method='ssm'
        while error>tol:
            if method=='ssm':
                R_new,K_new=self._successive_sub()
                frac_new=self.phases[0]['frac']
                if abs(frac_new-frac_old)<.1 and 1e-5<sum([(i-1)**2 for i in R_new])<1e-3 and 0<frac_new<1\
                        and sum([(i-1)**2 for i in R_new])/sum([(i-1)**2 for i in R_old])>.8:
                    method='assm'
                frac_old=frac_new
                R_old,K_old=R_new,K_new
            elif method=='assm':
                R_new,K_new=self._a_successive_sub(R_old,K_old)
                if sum([(i-1)**2 for i in R_new])>sum([(i-1)**2 for i in R_old]):
                    method='ssm'
                R_old=R_new
            else:
                assert False
            error=sum([(i-1)**2 for i in R_new])
        self._split_aq()
        return self.phases