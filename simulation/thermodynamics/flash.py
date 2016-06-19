import math
import peng_robinson
import numpy as np
from scipy.optimize import minimize

R=8.314
ln=math.log
# Pre-processing
class Prep:
    def __init__(self,pvt,P,T):
        self.P=P
        self.T=T
        self.pvt=self.set_types(pvt)
        self.tangents=self.set_binary_tangent_planes()
        self.dominant=self.set_dominant_pair()
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
            pvt_ij = {i: {'params': pvt[i]['params'], 'comp': comp_i, 'k':pvt[i]['k']}, j: {'params': pvt[j]['params'], 'comp': comp_j, 'k':pvt[j]['k']}}
            fug_ij= peng_robinson.fug_minimum_gibbs(pvt_ij,T, P, 1, phase='heavy')
            fug_i0 = peng_robinson.fug_minimum_gibbs({i:{'params': pvt[i]['params'], 'comp': 1., 'k':pvt[i]['k']}},T0,P0, 1, phase='heavy')[i]
            fug_j0 = peng_robinson.fug_minimum_gibbs({j:{'params': pvt[j]['params'], 'comp': 1., 'k':pvt[j]['k']}},T0,P0, 1, phase='heavy')[j]
            return  sign*(comp_i*(ln(comp_i)+fug_ij[i]-fug_i0)+comp_j*(ln(comp_j)+fug_ij[j]-fug_j0))
        def objective(x,P,T,i,j,get_ideal,sign):
            return func(x, P, T, i, j, sign)-get_ideal(x,sign)
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
    def set_dominant_pair(self):
        pvt=self.pvt
        tangents=self.tangents
        max_moles=0
        dominant=None
        for i in pvt:
            for j in pvt:
                if j==i:
                    pass
                elif len(tangents[(i,j)])==2 or tangents[(i,j)][0]['GRT']>0:
                    if pvt[i]['comp']+pvt[j]['comp']>max_moles:
                        dominant=(i,j)
                        max_moles=pvt[i]['comp']+pvt[j]['comp']
        return dominant
    def set_types(self,pvt):
        T=self.T
        P=self.P
        fug_h=peng_robinson.fug_minimum_gibbs(pvt,T,P,1.,phase='heavy')
        K_ideal={}
        comp_l={}
        for i in pvt:
            pvt_i=pvt[i]
            K_ideal_i=math.exp(peng_robinson.fug_minimum_gibbs({i:pvt_i},T,P,1.,comp={i:1.},phase='heavy')[i]
                     -peng_robinson.fug_minimum_gibbs({i:pvt_i},T,P,1.,comp={i:1.},phase='light')[i])
            K_ideal[i]=K_ideal_i
            comp_l[i]=K_ideal_i*pvt_i['comp']
            if K_ideal_i<1e-5:
                pvt[i]['type']="non-volatile"
        for i in pvt:
            if 'type' in pvt[i]:
                pass
            else:
                K=math.exp(fug_h[i]+ln(K_ideal[i])-peng_robinson.fug_minimum_gibbs(pvt,T,P,1.,comp=comp_l,phase='light')[i])
                if K>1e7:
                    pvt[i]['type']='non-condensible'
                else:
                    pvt[i]['type'] ='normal'
        return pvt
