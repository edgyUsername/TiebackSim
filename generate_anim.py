import simulation.thermodynamics.peng_robinson as pr
import simulation.thermodynamics.setup_comp_data as setup_comp_data
import math
ln=math.log
comps=('ethanol','water')
T=90
P=101.3e3
G_pure_l=(lambda pvt,T,P:{i:pr.fug_minimum_gibbs(setup_comp_data.get_props_and_equalize([(1e-6, j) if j!=i else (1,j) for j in pvt]),T,P,1.,phase='light') for i in pvt})(comps,T,P)
G_pure_h=(lambda pvt,T,P:{i:pr.fug_minimum_gibbs(setup_comp_data.get_props_and_equalize([(1e-6, j) if j!=i else (1,j) for j in pvt]),T,P,1.,phase='heavy') for i in pvt})(comps,T,P)

def G(pvt,T,P,phase):
    fugs = pr.fug_minimum_gibbs(pvt, T, P, 1., phase=phase)
    g = 0
    G_pure=G_pure_l if phase=='light' else G_pure_h
    for i in pvt:
        g += pvt[i]['comp'] * (G_pure[i][i] - fugs[i] - ln(pvt[i]['comp']))
    return g
x=[min(max(i/100.,1e-6),1-1e-6) for i in range(101)]
Gibbs_l=[G(setup_comp_data.get_props_and_equalize([(i, comps[0]),(1-i,comps[1])]),T,P,'light') for i in x]
Gibbs_h=[G(setup_comp_data.get_props_and_equalize([(i, comps[0]),(1-i,comps[1])]),T,P,'heavy') for i in x]
V=[i/100.for i in range(101)]
gl=''
gh=''
gt=''
H=[]
for i in range(101):
    gl+='(%s,%s,%s)'%(x[i], 1,Gibbs_l[i])
    gh += '(%s,%s,%s)' % (x[i], 0, Gibbs_h[i])
    for j in range(100):
        if not i==100 and i%10==0 and j%10==0:
            try:
                f=.8
                h=(f-x[i]*V[j])/(1-V[j])
                h2=(f-x[i+10]*V[j])
                h3 = (f - x[i+10] * V[j+10]) / (1 - V[j+10])
                h4 = (f- x[i] * V[j+10]) / (1 - V[j+10])
            except:
                pass
            else:
                if 0<h<1 and 0<h2<1 and 0<h3<1 and 0<h4<1:
                    gt += '(%s,%s,%s)' % (x[i], V[j], V[j]*Gibbs_l[i]+(1-V[j])*G(setup_comp_data.get_props_and_equalize([(h, comps[0]),(1-h,comps[1])]),T,P,'heavy'))
                    gt += '(%s,%s,%s)' % (x[i+10], V[j], V[j]*Gibbs_l[i+10]+(1-V[j])*G(setup_comp_data.get_props_and_equalize([(h2, comps[0]),(1-h,comps[1])]),T,P,'heavy'))
                    gt += '(%s,%s,%s)' % (x[i+10], V[j+10], V[j+10]*Gibbs_l[i+10]+(1-V[j+10])*G(setup_comp_data.get_props_and_equalize([(h3, comps[0]),(1-h,comps[1])]),T,P,'heavy'))
                    gt += '(%s,%s,%s)' % (x[i], V[j+10], V[j+10]*Gibbs_l[i]+(1-V[j+10])*G(setup_comp_data.get_props_and_equalize([(h4, comps[0]),(1-h,comps[1])]),T,P,'heavy'))

wrapper="""
\\begin{figure}
	\\centering
	\\begin{tikzpicture}
	\\begin{axis}[width=7.8cm,height=6.2cm,view={60}{30},
	title={Dimensionless Gibbs Plot},
	xlabel={$x_{\mathrm{ethanol}}$},
	ylabel={$V$},
	zlabel={$\Delta G/RT$},
	zmax=1,zmin=-1]
	\\addplot3[patch,colormap/cool,patch type=rectangle]
	coordinates{%s};"""%gt
b="""
    \\addplot3
	coordinates{%s};
	\\addplot3
	coordinates{%s};
	\\end{axis}\\end{tikzpicture}\\end{figure}
"""%(gh,gl)
print gt
print gh
print gl
pass