import p_and_a as pa
import copy
calc=pa.calcPressureDrop
mean={'m':8.,'P0': 1796712.5441856971, 'v_g': 1.2787357748465502, 'T': 4.080327374649706, 'v_l': 0.55583613436757029, 'm_l': 0.017974905736765361, 'pattern': 'elongated bubble', 'm_g': 1.080278388566883e-05, 'hold_up': 0.40088432954763276, 'enthalpy': -1759441.5196773682, 'r_l': 594.04027654902427, 'geo_index': 99900, 'Cp_l': 1825.3164389077942, 'Cp_g': 2109.3699874156773, 'vapQ': .01, 'r_g': 14.284462925222096,'z':.1,'roughness':0.0000254,'ID':0.292,'L':1}
var={'m_l':1.1,'m_g':1.05,'r_g':1.03,'r_l':1.05,'roughness':10.,'z':1.01,'vapQ':1/1.04}
meanP=calc(mean['P0'],mean['m'],mean['vapQ'],mean['r_g'],mean['r_l'],mean['m_g'],mean['m_l'],mean['z'],mean['roughness'],mean['ID'],mean['L'],mean['T'])[0]
full_var=(meanP-calc(mean['P0'],mean['m'],mean['vapQ']*(var['vapQ']),mean['r_g']*(var['r_g']),mean['r_l']*(var['r_l']),mean['m_g']*(var['m_g']),mean['m_l']*(var['m_l']),mean['z']*(var['z']),mean['roughness']*(var['roughness']),mean['ID'],mean['L'],mean['T'])[0])/meanP
data=[]
for i in var:
    dat=copy.deepcopy(mean)
    dat[i]=dat[i]*var[i]
    P = calc(dat['P0'], dat['m'], dat['vapQ'], dat['r_g'], dat['r_l'], dat['m_g'], dat['m_l'],
                         dat['z'], dat['roughness'], dat['ID'], dat['L'], dat['T'])[0]
    data.append((i,-(100*(P-meanP)/(meanP*full_var))))
s=sum([i[1] for i in data])
data=[(i[0],100*abs(i[1])/s) for i in data]
for i in data:
    print i
