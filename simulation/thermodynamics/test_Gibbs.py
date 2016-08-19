import setup_comp_data
import copy
import peng_robinson as pr
import math
dc=copy.deepcopy
T,P=98,1e5
ln=math.log


pvt=setup_comp_data.get_props_and_equalize([(1,'water'),(1,'octane'),(1,'methanol')])
G_pure={'methanol':pr.fug_minimum_gibbs(setup_comp_data.get_props_and_equalize([(1e-6, 'water'), (1e-6, 'octane'),(1,'methanol')]),T,P,1.,phase='light') ,'octane':pr.fug_minimum_gibbs(setup_comp_data.get_props_and_equalize([(1e-6, 'water'), (1, 'octane'),(1e-6,'methanol')]),T,P,1.,phase='light'),'water':pr.fug_minimum_gibbs(setup_comp_data.get_props_and_equalize([(1, 'water'), (1e-6, 'octane'),(1e-6,'methanol')]),T,P,1.,phase='light')}
AB=[-1,1,G_pure['water']['water']-G_pure['octane']['octane']]
AC=[-1,0,G_pure['methanol']['methanol']-G_pure['octane']['octane']]
matrix=[AB[1]*AC[2]-AB[2]*AC[1],AB[2]*AC[0]-AB[0]*AC[2],AB[0]*AC[1]-AB[1]*AC[0]]
d=-G_pure['water']['water']*matrix[2]-matrix[1]
matrix.append(d)
def G(pvt):
    fugs=pr.fug_minimum_gibbs(pvt,T,P,1.,phase='light')
    # g=0
    # for i in pvt:
    #     g+=fugs[i]+ln(pvt[i]['comp'])
    # return g
    g=0
    for i in pvt:
        g+=pvt[i]['comp']*(G_pure[i][i]-fugs[i]-ln(pvt[i]['comp']))
    return g
coords=[]
c3=[]
for i in range(11):
    oct=i/10.
    c1=[]
    c2=[]
    for j in range(11):
        wat=j/10.
        if oct+wat==1:
            if oct==1:
                pvt = setup_comp_data.get_props_and_equalize(
                    [(1e-6, 'water'), (oct, 'octane'), (1e-6, 'methanol')])
            elif oct==0:
                pvt = setup_comp_data.get_props_and_equalize(
                    [(wat, 'water'), (1e-6, 'octane'), (1e-6, 'methanol')])
            else:
                pvt = setup_comp_data.get_props_and_equalize(
                    [(wat, 'water'), (oct, 'octane'), (1e-6, 'methanol')])
        elif oct+wat<1:
            if oct+wat==0:
                pvt = setup_comp_data.get_props_and_equalize(
                    [(1e-6, 'water'), (1e-6, 'octane'), (1 - oct - wat, 'methanol')])
            elif oct==0:
                pvt = setup_comp_data.get_props_and_equalize(
                    [(wat, 'water'), (1e-6, 'octane'), (1 - oct - wat, 'methanol')])
            elif wat==0:
                pvt = setup_comp_data.get_props_and_equalize([(1e-6, 'water'), (oct, 'octane'), (1 - oct - wat, 'methanol')])
            else:
                pvt = setup_comp_data.get_props_and_equalize([(wat, 'water'), (oct, 'octane'),(1-oct-wat,'methanol')])
                # c1.append((oct, wat, G(pvt)))
        if oct+wat<=1:
            fugsl,fugsh=pr.fug_minimum_gibbs(pvt,T,P,1.,phase='light'),pr.fug_minimum_gibbs(pvt,T,P,1.,phase='heavy')
            c2.append((oct,wat,1-oct-wat,pr.find_vapor_frcn(pvt,{z:fugsh[z]-fugsl[z] for z in fugsl})[0]))
            c1.append((oct,wat,G(pvt)))
    coords.append(c1)
    c3.append(c2)
cmin=(0,0)
count0=0
for i in coords:
    count1=0
    for j in i:
        if j[2]<coords[cmin[0]][cmin[1]][2]:
            cmin=(count0,count1)
        count1+=1
    count0+=1
print coords[cmin[0]][cmin[1]],'\n'
c=''
#-oct*G_pure['octane']['octane']-wat*G_pure['water']['water']-(1-oct-wat)*G_pure['methanol']['methanol']
for i in range(len(coords)-1):
    for j in range(len(coords)-1):
        if j+1<len(coords[i+1]) and j+1<len(coords[i]):
            c+=str(coords[i][j])+str(coords[i][j+1])+str(coords[i+1][j+1])+str(coords[i+1][j])
print c
d=''
for i in range(len(c3)-1):
    for j in range(len(c3)-1):
        if j+1<len(c3[i+1]) and j+1<len(c3[i]):
            d+=str(c3[i][j])+str(c3[i][j+1])+str(c3[i+1][j+1])+str(c3[i+1][j])
print d
pass
