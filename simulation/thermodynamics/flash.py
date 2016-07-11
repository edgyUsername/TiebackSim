import math
import peng_robinson
import numpy as np
from scipy.optimize import minimize,minimize_scalar

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
def normalize_comps(comp):
    total_moles = sum(comp.values())
    for i in comp:
        comp[i]= comp[i]/total_moles
    return comp
# Pre-processing
class Prep:
    def __init__(self,pvt,P,T):
        self.P=P
        self.T=T
        self.pvt=self.set_types(pvt.copy())
        self.tangents=self.set_binary_tangent_planes()
    def set_binary_tangent_planes(self):
        """ binary tangents assuming a liquid phase
        :return: dictionary of binary tangent planes
        """
        pvt= self.pvt
        P= self.P
        T= self.T
        def ideal(P,T,i,j):
            f_start=func([.01,.99],P,T,i,j,1)
            f_end=func([.99,.01],P,T,i,j,1)
            grad=(f_end-f_start)/.98
            intercept=f_start-grad*.01
            def find_ideal(x,sign):
                return sign*(intercept+x[0]*grad)
            return find_ideal
        def func(x,P,T,i,j,sign):
            comp_i,comp_j=x
            P0, T0 = 1e5, 25
            pvt_ij = {i: {'params': pvt[i]['params'].copy(), 'comp': comp_i, 'k':pvt[i]['k']}, j: {'params': pvt[j]['params'].copy(), 'comp': comp_j, 'k':pvt[j]['k']}}
            fug_ij= peng_robinson.fug_minimum_gibbs(pvt_ij,T, P, 1, phase='heavy')
            fug_i0 = peng_robinson.fug_minimum_gibbs({i:{'params': pvt[i]['params'].copy(), 'comp': 1., 'k':pvt[i]['k']}},T0,P0, 1, phase='heavy')[i]
            fug_j0 = peng_robinson.fug_minimum_gibbs({j:{'params': pvt[j]['params'].copy(), 'comp': 1., 'k':pvt[j]['k']}},T0,P0, 1, phase='heavy')[j]
            return  sign*(comp_i*(ln(comp_i)+fug_ij[i]-fug_i0)+comp_j*(ln(comp_j)+fug_ij[j]-fug_j0))
        def objective(x,P,T,i,j,get_ideal,sign):
            return func(x, P, T, i, j, sign)-get_ideal(x,sign)
        for i in pvt:
            for j in pvt:
                if not j==i:
                    print 'x_{%s}&\\frac{\Delta\\bar{G}_{\mathrm{mix}}}{RT}'%i
                    x=''
                    for n in range(1,100):
                        c=n/100.
                        ans=objective([c,1-c],P,T,i,j,ideal(P,T,i,j),1)
                        x+=str((c,ans))
                    print x
        cons = ({'type': 'eq',
                 'fun': lambda x: np.array([1-x[0]-x[1]]),
                 'jac': lambda x: np.array([-1.,-1.])},
                {'type': 'ineq',
                 'fun': lambda x: np.array([x[0]-1e-3]),
                 'jac': lambda x: np.array([1., 0.])})
        tangents={}
        delta=1e-2
        tol=1e-2
        for i in pvt:
            for j in pvt:
                get_ideal=ideal(P,T,i,j)
                if not j==i and not (i,j) in tangents:
                    min_GRT = []
                    min_x=[]
                    sign=1
                    bounds=([delta/2,1-delta/2],[delta/2,1-delta/2])
                    c,error=0,tol+1
                    init1,init2=[delta, 1-delta],[1-delta, delta]
                    while error>tol and c<10:
                        val1 = minimize(objective, init1, args=(P,T,i,j,get_ideal,sign),
                                       constraints = cons, method = 'SLSQP',bounds=bounds)
                        val2 = minimize(objective, init2, args=(P,T,i,j,get_ideal,sign),
                                        constraints= cons, method= 'SLSQP',bounds=bounds)
                        error=abs(val1.x[0]-val2.x[0])
                        if error<tol and sign==1:
                            min_GRT.append(val1.fun)
                            min_x.append(val1.x)
                        elif sign==1:
                            min_GRT.append(val1.fun)
                            min_x.append(val1.x)
                            min_GRT.append(val2.fun)
                            min_x.append(val2.x)
                        c+=1
                        sign*=-1
                        bounds=([val1.x[0],val2.x[0]],[1-val2.x[0],1-val1.x[0]])
                        if val1.x[0]+delta>=val2.x[0]:
                            break
                        init1,init2=[val1.x[0]+delta,1-val1.x[0]-delta],[val2.x[0]-delta,1-val2.x[0]+delta]
                    c=0
                    tan_x=[]
                    tan_x2=[]
                    tan_GRT=[]
                    min_glob_1=None
                    min_glob_2=None
                    for point in min_GRT:
                        if c==0:
                            tan_GRT.append(point)
                            tan_x.append(min_x[c][0])
                            tan_x2.append(min_x[c][1])
                            min_glob_1=point
                        elif c==1:
                            tan_GRT.append(point)
                            tan_x.append(min_x[c][0])
                            tan_x2.append(min_x[c][1])
                            min_glob_2=point
                        elif point<min_glob_1:
                            tan_GRT[0]=point
                            tan_x[0]=min_x[c][0]
                            tan_x2[0]=min_x[c][1]
                            min_glob_1 = point
                        elif point<min_glob_2:
                            tan_GRT[1] = point
                            tan_x[1] = min_x[c][0]
                            tan_x2[1] = min_x[c][1]
                            min_glob_2 = point
                        c+=1
                    tangents[(i,j)]=[{'x':tan_x[z],'GRT':tan_GRT[z]} for z in range(len(tan_GRT))]
                    tangents[(j,i)]=[{'x':tan_x2[z],'GRT':tan_GRT[z]} for z in range(len(tan_GRT))]
        return tangents

    def set_types(self,pvt):
        T=self.T
        P=self.P
        fug_h=peng_robinson.fug_minimum_gibbs(pvt,T,P,1.,phase='heavy')
        K_ideal={}
        comp_l={}
        for i in pvt:
            pvt_i=pvt[i].copy()
            K_ideal_i=exp(peng_robinson.fug_minimum_gibbs({i:pvt_i},T,P,1.,comp={i:1.},phase='heavy')[i]
                     -peng_robinson.fug_minimum_gibbs({i:pvt_i},T,P,1.,comp={i:1.},phase='light')[i])
            K_ideal[i]=K_ideal_i
            comp_l[i]=K_ideal_i*pvt_i['comp']
            if K_ideal_i<1e-5:
                pvt[i]['type']="non-volatile"
        for i in pvt:
            if 'type' in pvt[i]:
                pass
            else:
                K=exp(fug_h[i]+ln(K_ideal[i])-peng_robinson.fug_minimum_gibbs(pvt,T,P,1.,comp=comp_l,phase='light')[i])
                if K>1e7:
                    pvt[i]['type']='non-condensible'
                else:
                    pvt[i]['type'] ='normal'
        return pvt

class Phase:
    def __init__(self,pvt,P,T,tangents,mole_fraction):
        self.pvt=pvt.copy()
        self.T=T
        self.P=P
        self.tangents=tangents
        self.dominant=self.set_dominant_pair()
        self.moles=mole_fraction

    def set_dominant_pair(self):
        pvt = self.pvt
        tangents = self.tangents
        max_moles = 0
        dominant = None
        for i in pvt:
            for j in pvt:
                if j == i:
                    pass
                elif len(tangents[(i, j)]) == 2 or tangents[(i, j)][0]['GRT'] > 0:
                    if pvt[i]['comp'] + pvt[j]['comp'] > max_moles:
                        dominant = (i, j)
                        max_moles = pvt[i]['comp'] + pvt[j]['comp']
        return dominant
    def get_rel_solubilities(self):
        pvt=self.pvt
        T=self.T
        P=self.P
        dominant=self.dominant
        tangents=self.tangents
        r_1={}
        r_2={}
        for i in pvt:
            if i in dominant:
                pass
            else:
                miscible=len(tangents[(i,dominant[0])])==1 and len(tangents[(i,dominant[1])])==1
                if miscible:
                    fug_i=peng_robinson.fug_minimum_gibbs({i:pvt[i].copy()},T,P,1.,comp={i:1.},phase="heavy")[i] #ln(fugacity coef of i)
                    fug_1=peng_robinson.fug_minimum_gibbs({dominant[0]:pvt[dominant[0]].copy()},T,P,1.,comp={dominant[0]:1.},phase="heavy")[dominant[0]]
                    fug_2=peng_robinson.fug_minimum_gibbs({dominant[1]:pvt[dominant[1]].copy()},T,P,1.,comp={dominant[1]:1.},phase="heavy")[dominant[1]]
                    si_1=fug_i+ ln(pvt[i]['comp']) - ln(1-pvt[i]['comp'])-fug_1
                    si_2=fug_i+ ln(pvt[i]['comp']) - ln(1-pvt[i]['comp'])-fug_2

                    r_1[i]=si_1/(si_1+si_2)
                    r_2[i]=si_2/(si_1+si_2)
                else:
                    D1=abs(tangents[(i,dominant[0])][0]['x']-tangents[(i,dominant[0])][1]['x']) if len(tangents[(i,dominant[0])])==2 else 0
                    D2=abs(tangents[(i,dominant[1])][0]['x']-tangents[(i,dominant[1])][1]['x']) if len(tangents[(i,dominant[1])])==2 else 0

                    r_1[i] = (1-D1)/((1-D1)+(1-D2)/(1-D1))
                    r_2[i] = (1-D2)/((1-D2)+(1-D1)/(1-D2))

        return (r_1,r_2)
    def split_init(self):
        pvt=self.pvt
        dominant=self.dominant
        if dominant==None:
            comp={}
            for i in pvt:
                comp[i]=pvt[i]['comp']
            return [comp]
        else:
            dom_tangents=self.tangents[dominant]
            comp_1,comp_2={},{}
            r_1,r_2= self.get_rel_solubilities()
            for i in pvt:
                if not i in dominant:
                    comp_1[i]=r_1[i]*pvt[i]['comp']
                    comp_2[i]=r_2[i]*pvt[i]['comp']
                elif i==dominant[0]:
                    comp_1[i]=max(dom_tangents[0]['x'],dom_tangents[1]['x'])
                    comp_2[i]=min(dom_tangents[0]['x'],dom_tangents[1]['x'])
                elif i==dominant[1]:
                    comp_1[i] = 1-max(dom_tangents[0]['x'], dom_tangents[1]['x'])
                    comp_2[i] = 1-min(dom_tangents[0]['x'], dom_tangents[1]['x'])
                else:
                    assert  False, 'component "%s" not found in system' %i
            comp_1,comp_2=normalize_comps(comp_1),normalize_comps(comp_2)
            return (comp_1,comp_2)
    def split_phase(self):
        pvt=self.pvt
        comps=self.split_init()
        if len(comps)==1:
            new_pvt={}
            for i in comps[0]:
                new_pvt[i]=pvt[i].copy()
                new_pvt[i]['comp']=comps[0][i]
            return [Phase(new_pvt,self.P,self.T,self.tangents,self.moles)]
        else:
            index=[]
            for i in pvt:
                index.append(i)
            comp_1=[comps[0][j] for j in index]
            comp_2=[comps[1][j] for j in index]
            z=[pvt[j]['comp'] for j in index]
            def func(frac):
                return ((1-frac)**2)*np.dot(comp_2,comp_2)+2*frac*(1-frac)*np.dot(comp_2,comp_1)+(frac**2)*np.dot(comp_1,comp_1)-np.dot(z,z)\
                       -2*(1-frac)*np.dot(z,comp_2)-2*frac*np.dot(z,comp_1)
            bounds=[0.,1.]
            fraction=minimize_scalar(func,bounds=bounds,method='Bounded').x
            new_pvt1={}
            new_pvt2={}
            c=0
            for i in index:
                new_pvt1[i],new_pvt2[i]=pvt[i].copy(),pvt[i].copy()
                new_pvt1[i]['comp'],new_pvt2[i]['comp']=comp_1[c],comp_2[c]
                c+=1
            return [Phase(new_pvt1,self.P,self.T,self.tangents,self.moles*fraction),
                    Phase(new_pvt2,self.P,self.T,self.tangents,self.moles*(1-fraction))]

class System:
    def __init__(self,pvt,P,T):
        self.P=P
        self.T=T
        prep=Prep(pvt,P,T)
        self.pvt=prep.pvt
        self.tangents=prep.tangents
        print pvt
        self._set_g_pure()
        self.stable=False
        self.split_liquid_phases=None
        self.liquid_phases=[Phase(self.pvt.copy(),P,T,self.tangents,1.)]
        self.phases=[]
    def _check_stability(self):
        def _make_new_phases(self):
            new_phases = []
            updated = False
            pos=0
            for phase in self.liquid_phases:
                sub_phases = phase.split_phase() if not updated else [phase]
                if len(sub_phases) == 1:
                    new_phases.append(sub_phases[0])
                else:
                    updated = True
                    new_phases.append(sub_phases[0])
                    new_phases.append(sub_phases[1])
                pos+=1

            self.split_liquid_phases = new_phases
        _make_new_phases(self)
        stable=len(self.liquid_phases)==len(self.split_liquid_phases)
        return stable
    def _set_g_pure(self):
        g_pure_i = {'heavy': {}, 'light': {}}
        T = self.T + 273.15
        P = self.P
        pvt = self.pvt
        for i in pvt:
            z_i_roots, eos_params = peng_robinson.solve_PR_for_Z(generate_pvt(pvt, [i]), (T - 273.15), P,comp={i: 1.})
            v_heavy = z_i_roots[0] * R * T / P
            v_light = z_i_roots[-1] * R * T / P
            a, b = eos_params['comp_params'][i]['a_i'], eos_params['comp_params'][i]['b_i']
            g_pure_i['heavy'][i]=P * v_heavy / (R * T) + ln(R * T / (v_heavy - b)) + (a / (R * T * b * (8 ** .5)))* ln(
                (2 * v_heavy + 2 * b - b * (8 ** .5)) / (2 * v_heavy + 2 * b + b * (8 ** .5)))
            g_pure_i['light'][i]=P * v_light / (R * T) + ln(R * T / (v_light - b)) + (a / (R * T * b * (8 ** .5)))* ln(
                (2 * v_light + 2 * b - b * (8 ** .5)) / (2 * v_light + 2 * b + b * (8 ** .5)))
        self.g_pure=g_pure_i
    def _g_peng_robinson(self, comp, phase='heavy'):
        T, P, pvt,g_pure= self.T, self.P, self.pvt,self.g_pure
        index=[i for i in comp]
        n=len(index)
        Z_roots, eos_params = peng_robinson.solve_PR_for_Z(pvt, T, P, comp={index[i]: comp[index[i]] for i in range(n)})
        v = Z_roots[-1] * R * (T + 273.15) / P if phase == 'light' else Z_roots[0] * R * (T + 273.15) / P
        a, b = eos_params['a_m'], eos_params['b_m']
        T = T + 273.15
        return P * v / (R * T) + ln(R * T / (v - b)) + (a / (R * T * b * (8 ** .5))) * ln(
            (2 * 2 * b - b * (8 ** .5)) / (2 * v + 2 * b + b * (8 ** .5))) \
               + sum([comp[index[j]] * ln(comp[index[j]]) for j in range(n)]) - sum([comp[index[j]] * g_pure[phase][index[j]] for j in range(n)])

    def _flash_liquids(self,vapor=False):
        """updates phases in system  by checking current phases for stability and splitting unstable liquid phases"""
        if not vapor:
            stable=self._check_stability()
        else:
            stable=False
        if stable:
            self.stable=True
        else: # continues to minimize gibbs if not all stable liquid phases are found

            pvt=self.pvt
            P=self.P
            T=self.T
            index=[]
            gas=-1
            for i in pvt:
                index.append(i)
            if not vapor:
                x0=[l.moles for l in self.split_liquid_phases] #l1,l2,...,lk,x1l1,x2l1...xnl1...xnlk
                k=len(self.split_liquid_phases) # number of phases
                n=len(index) # number of components
                for l in self.split_liquid_phases:
                    pvt_l=l.pvt
                    for i in index:
                        x0.append(pvt_l[i]['comp'])
            else:
                x0=[l['phase'].moles for l in self.phases]
                k = len(self.phases)  # number of phases
                n = len(index)  # number of components
                c=0
                for l in self.phases:
                    if l['type']=='gas':
                        gas=c
                    pvt_l = l['phase'].pvt
                    for i in index:
                        x0.append(pvt_l[i]['comp'])
                    c+=1

            ########## find gibbs of mixing for each phase
            def _total_gibbs(x):
                x=[c if c>0 else 1e-19 for c in x]
                G_total=0
                for p in range(k):
                    frac_p=x[p]
                    comp_p={index[i]:x[k+n*p+i] for i in range(n)}
                    if vapor and p==gas:
                        fugij=peng_robinson.fug_minimum_gibbs(pvt,T,P,1,comp=comp_p,phase='light')
                        # G_total+=frac_p*self._g_peng_robinson(comp_p,phase='light')
                        G_total+=frac_p*sum([comp_p[i]*(ln(comp_p[i])+fugij[i]) for i in index])
                    else:
                        fugij = peng_robinson.fug_minimum_gibbs(pvt, T, P, 1, comp=comp_p, phase='heavy')
                        # G_total+=frac_p*self._g_peng_robinson(comp_p)
                        G_total += frac_p * sum([comp_p[i] * (ln(comp_p[i]) + fugij[i]) for i in index])
                return G_total
            def fug_diff(x):
                x = [c if c > 0 else 1e-18 for c in x]
                fug_eq_list=[]
                fug_p=[]
                comp_p=[]
                for p in range(k):
                    comp_p.append({index[j]:x[k+n*p+j] for j in range(n)})
                    if vapor and p==k-1:
                        fug_p.append(peng_robinson.fug_minimum_gibbs(pvt,T,P,1.,comp=comp_p[p],phase='light'))
                    else:
                        fug_p.append(peng_robinson.fug_minimum_gibbs(pvt, T, P, 1., comp=comp_p[p], phase='heavy'))
                for i in range(n):
                    fug_i = []
                    fug_diff_list = []
                    for p in range(k):
                        fug_i.append(ln(comp_p[p][index[i]])+fug_p[p][index[i]])
                    for p in range(k-1):
                        fug_diff_list.append(fug_i[p]-fug_i[p+1])
                    fug_eq_list+=fug_diff_list
                return fug_eq_list

            cons = ({'type': 'eq',
                     'fun': lambda x: np.array([1 - sum([x[l] for l in range(k)])]+[1-sum([x[k+i+n*p] for i in range(n)]) for p in range(k)]),
                     'jac': lambda x: np.array([[(-1. if c in range(k) else 0)for c in range(len(x))]]
                                               +[[(-1. if c in range(k+n*p,k+n*(p+1)) else 0) for c in range(len(x))]for p in range(k)])},
                    {'type': 'ineq',
                     'fun': lambda x: np.array([c for c in x]+[1-c for c in x]),
                     'jac': lambda x: np.array([[1. if c1==c2 else 0 for c1 in range(len(x))] for c2 in range(len(x))]
                                               +[[-1. if c1==c2 else 0 for c1 in range(len(x))] for c2 in range(len(x))])},
                    {'type':'eq',
                     'fun': lambda x: np.array([self.pvt[index[i]]['comp']-sum([x[k+i+n*p] for p in range(k)]) for i in range(n)])})


            # bounds=[[0.,1.] for c in x0]
            res=minimize(_total_gibbs,x0,method='SLSQP',constraints=cons)
            # parse results from minimizer
            new_phases=[]
            for p in range(k):
                fraction=res.x[p]
                pvt_k=generate_pvt(pvt,[i for i in index])
                for i in range(n):
                    pvt_k[index[i]]['comp']=res.x[(k-1)+(i+1)+n*p]
                new_phases.append(Phase(pvt_k,P,T,self.tangents,fraction))
            if not vapor:
                self.split_liquid_phases=new_phases
            else:
                self.phases=[{'type':'gas' if p==k-1 else 'liquid','phase':new_phases[p]} for p in range(k)]

    def _get_bubble_props(self):
        T=self.T
        P=self.P
        print T
        pvt=self.pvt.copy()
        liquids=self.liquid_phases
        volatiles=[]
        comps=[]
        index=[]
        k=0
        n=0
        for i in pvt:
            # remove non-volatiles
            if pvt[i]['type']!="non-volatile":
                volatiles.append(i)
                index.append(i)
                n+=1
        if len(volatiles)==0:
            return None
        else:
            pvt=generate_pvt(pvt,volatiles)
            for l in liquids:
                k+=1
                new_pvt=generate_pvt(l.pvt,volatiles)
                for i in index:
                    comps.append(new_pvt[i]['comp'])
            def _get_lnK(T,comps_y=None):
                comp_x={}
                lnK=[]
                for i in range(len(index)):
                    comp_x[index[i]] = pvt[index[i]]['comp']
                if not comps_y==None:
                    comp_y={}
                    for i in range(len(comps_y)):
                        comp_y[index[i]] = comps_y[i] if comps_y[i]>1e-18 else 1e-18
                    fug_l=peng_robinson.fug_minimum_gibbs(pvt,T,P,1.,comp=comp_x,phase='heavy')
                    fug_v=peng_robinson.fug_minimum_gibbs(pvt,T,P,1.,comp=comp_y,phase='light')

                    for i in index:
                        lnK.append(fug_l[i]-fug_v[i])
                else:
                    for i in index:
                        pvt_i = pvt[i].copy()
                        lnK.append(peng_robinson.fug_minimum_gibbs({i: pvt_i}, T, P, 1., comp={i: 1.}, phase='heavy')[i]
                            - peng_robinson.fug_minimum_gibbs({i: pvt_i}, T, P, 1., comp={i: 1.}, phase='light')[i])
                return lnK
            def func(x):
                lnK=[_get_lnK(x[-1],comps_y=x[:-1]) for p in range(k)]
                return sum([sum([(exp(lnK[p][i])*comps[i+n*p]-x[i])**2 for i in range(n)]) for p in range(k)])+(1-sum([x[i] for i in range(n)]))**2
            # initial estimates
            x0=[]
            lnK_init=_get_lnK(T)
            tot=0
            for i in range(len(index)):
                y_i=math.exp(lnK_init[i])*pvt[index[i]]['comp']
                tot+=y_i
                x0.append(y_i)
            x0=[i/tot for i in x0]
            # for t in range(-273,130):
            #     func(x0+[t])
            # initial temp estimate
            T01,T02=T,T
            f1, f2 = func(x0 + [T01]), func(x0 + [T02])
            while (f1<.1 or (f1>.9 and f1<1.2)) and (f2<.1 or (f2>.9 and f2<1.2)):
                T01+=10
                T02=T02-20 if T02>-253 else T02*.2
                f1,f2=func(x0+[T01]),func(x0+[T02])
                if T02<-253 and (f2<.1 or f2>.9):
                    return None
            T0=T01 if f1>=.1 else T02
            x0.append(T0)
            print x0
            bounds=[[0.,1.]for i in index]+[[-273,None]]
            res=minimize(func,x0,method='SLSQP',bounds=bounds)
            return {'comp':normalize_comps({index[i]:res.x[i] for i in range(len(index))}),'T':res.x[-1]}
    def _equilibriate(self):

        while not self.stable:
            self._flash_liquids()
            new_phases=[]
            for p in self.split_liquid_phases:
                collapse=p.moles<1e-3
                if not collapse:
                    new_phases.append(p)
                else:
                    self.stable=True
            if not self.stable:
                self.liquid_phases=new_phases
            self.split_liquid_phases=None
