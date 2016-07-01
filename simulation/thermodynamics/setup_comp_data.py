import csv
import sys
import peng_robinson
import math
from scipy.optimize import minimize
from os import path
sys.path.append( path.dirname( path.dirname(path.dirname( path.abspath(__file__)) ) ) )
exp=math.exp
ln=math.log
sinh=math.sinh
cosh=math.cosh
def all_crit_properties():
    """returns dict of all properties from Perrys_crit.csv"""
    properties={}
    with open(sys.path[-1]+'/databases/Perrys_crit.csv','rb') as csvfile:
        prop_reader=csv.DictReader(csvfile)
        for row in prop_reader:
            if row['mm'] =="":
                pass
            else:
                names=[row['name'].replace(' ','').replace(',','').replace('-','').lower(),]
                for i in row['alt_name'].split(','):
                    names.append(i.replace(' ','').replace('-','').lower()) if not i.replace(' ','').replace('-','') =="" else None
                for n in names:
                    properties[n]={
                        'p_c':float(row['p_c']),				                                    #bara
                        't_c':float(row['t_c']),				                                    #K
                        'v_c':float(row['v_c']),				                                    #m3/kmol
                        'acc':float(row['acc']),				                                    #accentric factor
                        'mm':float(row['mm']),					                                    #kg/mol
                        '_id':row['name'].replace(' ','').replace(',','').replace('-','').lower(),
                        'cl':[float(row['c%s'%(i+1)]) for i in range(5)],                           #J/(kmol K)
                        'cg':[float(row['cg%s'%(i+1)]) for i in range(5)],                          #J/(kmol K)
                        'hf':float(row['hf'])                                                       #J/mol
                    }
                    if row['eq']=='1':
                        properties[n]['heat_l']={'C':lambda c,T,t_c: sum([c[i]*(T+273.15)**i for i in range(5)]),
                                               'h':lambda c,T,t_c: (T-25)*(sum([c[i]*(T+273.15)**i for i in range(5)])+sum([c[i]*(25+273.15)**i for i in range(5)]))/2} # sum([c[i]*((T+273.15)**(i+1))/(i+1) for i in range(5)])-sum([c[i]*((15+273.15)**(i+1))/(i+1) for i in range(5)])}
                    elif row['eq']=='2':
                        properties[n]['heat_l']={'C':lambda c,T,t_c: sum([-(-c[4]**2/5)*((T+273.15)**5)/(t_c**5),
                                                         (5*(-c[4]**2/5)+(-c[2]*c[3]/2))*((T+273.15)**4)/(t_c**4),
                                                         (-10*(-c[4]**2/5)-4*(-c[2]*c[3]/2)-(-c[2]**2/3))*((T+273.15)**3)/(t_c**3),
                                                         (10*(-c[4]**2/5)+6*(-c[2]*c[3]/2)+3*(-c[2]**2/3)+(-c[0]*c[3]))*((T+273.15)**2)/(t_c**2),
                                                         (-5*(-c[4]**2/5)-4*(-c[2]*c[3]/2)-3*(-c[2]**2/3)-2*(-c[0]*c[3])-(-2*c[0]*c[2]))*((T+273.15))/t_c,
                                                         (5+(c[1])),
                                                         (c[0]**2)*t_c/(t_c-(T+273.15))]),
                                               'h':lambda c,T,t_c: (T-25)*(sum([-(-c[4]**2/5)*((T+273.15)**5)/(t_c**5),
                                                         (5*(-c[4]**2/5)+(-c[2]*c[3]/2))*((T+273.15)**4)/(t_c**4),
                                                         (-10*(-c[4]**2/5)-4*(-c[2]*c[3]/2)-(-c[2]**2/3))*((T+273.15)**3)/(t_c**3),
                                                         (10*(-c[4]**2/5)+6*(-c[2]*c[3]/2)+3*(-c[2]**2/3)+(-c[0]*c[3]))*((T+273.15)**2)/(t_c**2),
                                                         (-5*(-c[4]**2/5)-4*(-c[2]*c[3]/2)-3*(-c[2]**2/3)-2*(-c[0]*c[3])-(-2*c[0]*c[2]))*((T+273.15))/t_c,
                                                         (5+(c[1])),
                                                         (c[0]**2)*t_c/(t_c-(T+273.15))])+sum([-(-c[4]**2/5)*((25+273.15)**5)/(t_c**5),
                                                         (5*(-c[4]**2/5)+(-c[2]*c[3]/2))*((25+273.15)**4)/(t_c**4),
                                                         (-10*(-c[4]**2/5)-4*(-c[2]*c[3]/2)-(-c[2]**2/3))*((25+273.15)**3)/(t_c**3),
                                                         (10*(-c[4]**2/5)+6*(-c[2]*c[3]/2)+3*(-c[2]**2/3)+(-c[0]*c[3]))*((25+273.15)**2)/(t_c**2),
                                                         (-5*(-c[4]**2/5)-4*(-c[2]*c[3]/2)-3*(-c[2]**2/3)-2*(-c[0]*c[3])-(-2*c[0]*c[2]))*((25+273.15))/t_c,
                                                         (5+(c[1])),
                                                         (c[0]**2)*t_c/(t_c-(25+273.15))]))/2}
                    if row['eqg']=='1':
                        properties[n]['heat_g']={'C':lambda c,T: sum([c[0],
                                                                     c[1]*(c[2]/((T+273.15)*sinh(c[2]/(T+273.15))))**2,
                                                                     c[3]*(c[4]/((T+273.15)*cosh(c[4]/(T+273.15))))**2]),
                                                 'h':lambda c,T: (T-25)*sum([sum([c[0],
                                                                                 c[1]*(c[2]/((T+273.15)*sinh(c[2]/(T+273.15))))**2,
                                                                                 c[3]*(c[4]/((T+273.15)*cosh(c[4]/(T+273.15))))**2]),
                                                                             sum([c[0],
                                                                                 c[1]*(c[2]/((25+273.15)*sinh(c[2]/(25+273.15))))**2,
                                                                                 c[3]*(c[4]/((25+273.15)*cosh(c[4]/(25+273.15))))**2])])/2}
                    elif row['eqg'] == '2':
                        properties[n]['heat_g'] = {'C': lambda c,T: sum([c[i]*(T+273.15)**i for i in range(5)]),
                                                   'h': lambda c,T: (T-25)*(sum([c[i]*(T+273.15)**i for i in range(5)])+sum([c[i]*(25+273.15)**i for i in range(5)]))/2}
    return properties

def get_props_and_equalize(compData):
    """compData is list with tuples (feed comp, comp name)"""
    ref=all_crit_properties()
    total_moles=0
    for i in compData:
        total_moles+=float(i[0])
    pvtParams={}
    for i in compData:
        try:
            compParams=ref[i[1].replace(',','').replace('-','').replace(' ','').lower()]
        except:
            print '"%s" not available in library or missing some params. You can ammend the Perrys_crit.csv file in the "databases" folder' %i[1]
            raise
        if compParams['_id'] in pvtParams:
            assert False, "components must be unique!"
        pvtParams[compParams['_id']]={'params':compParams,'comp':float(i[0])/total_moles}
    for i in pvtParams:
        pvtParams[i]['k'] = {j: 0 for j in pvtParams}
        pvtParams[i]['params']['t_b']=get_bp(pvtParams[i])
    return pvtParams

def get_bp(pvt_i):
    T = 25.
    P = 1.013e5
    i=pvt_i['params']['_id']
    pvt = {i:pvt_i}
    def _get_lnK(T, comps_y=None):
        comp_x = {i:1.}
        if not comps_y == None:
            comp_y = {i:comps_y[0]}
            fug_l = peng_robinson.fug_minimum_gibbs(pvt, T, P, 1., comp=comp_x, phase='heavy')[i]
            fug_v = peng_robinson.fug_minimum_gibbs(pvt, T, P, 1., comp=comp_y, phase='light')[i]
            lnK=fug_l - fug_v
        else:
            pvt_i = pvt[1].copy()
            lnK=peng_robinson.fug_minimum_gibbs({1: pvt_i}, T, P, 1., comp={i: 1.}, phase='heavy')[i]\
                - peng_robinson.fug_minimum_gibbs({1: pvt_i}, T, P, 1., comp={i: 1.}, phase='light')[i]
        return lnK

    def func(x):
        lnK = _get_lnK(x[-1], comps_y=x[:-1])
        return (exp(lnK) - x[0]) ** 2 + (1 - x[0]) ** 2

    x0 = [1.]
    # initial temp estimate
    T01, T02 = T, T
    f1, f2 = func(x0 + [T01]), func(x0 + [T02])
    while (f1 < .1 or (f1 > .9 and f1 < 1.2)) and (f2 < .1 or (f2 > .9 and f2 < 1.2)):
        T01 += 10
        T02 = T02 - 20 if T02 > -253 else T02 * .2
        f1, f2 = func(x0 + [T01]), func(x0 + [T02])
        if T02 < -253 and (f2 < .1 or f2 > .9):
            return None
    T0 = T01 if f1 >= .1 else T02
    x0.append(T0)
    bounds = [[0., 1.],[-273, None]]
    res = minimize(func, x0, method='SLSQP', bounds=bounds)
    return res.x[-1]
