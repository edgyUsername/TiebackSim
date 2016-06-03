import sys
from os import path
sys.path.append(path.dirname( path.abspath(__file__)) ) 
import settings
from simulation.fluid_flow import b_and_b as bb
from simulation.fluid_flow import p_and_a as pa
from simulation import simulator as sim


sysParams= {
    'P0':5515805.6,   #Pa abs
    'T0':79.44,        #deg C
    'Ql':317.9746,      #sm3/day
    'sysVars':[
        {
        'type':'Pipe',  #Pipe or Leg
        'params':(
            5000,       #length (m) 
            0.343,      #diameter (m)
            0.0000254,  #roughness (m)
            0.0127,     #thickness (m)
            0.0,        #inclination (deg)
            10.0        #U (W/m2/K)
            )
        },{
        'type':'Pipe',  
        'params':(
            5000,       
            0.343,      
            0.0000254,  
            0.0127, 
            -0.1718875964,
            10
            )
        },{
        'type':'Pipe',  
        'params':(
            200,       
            0.343,      
            0.0000254,  
            0.0127, 
            90,
            10
            )
        }
    ],
    'PVT':{								#black oil pvt data at wellhead
    'r_g':41.6468,						#density of gas (kg/m3)
    'r_l':799.2982,						#density of liquid
    'm_g':0.0000131,					#viscosity of gas (Pa s)
    'm_l':0.002,						#viscosity of liquid
    'vapQ':	0.8226974806477881			#mass vapour quality 
    }
}
"""
1bpd=1.840130787037037e-06		                	m3/s 
1scfd=3.2774128000000003e-07		                m3/s
1 lb/ft3=16.018					                  	kg/m3
1"=0.0254			    		                  	m
1psi=6894.757						                Pa
"""

sys=settings.newSystem(sysParams)
#sys.printGeometry()
#sim.simulate(sys,mode='devision',steps=300,devMultiplier=5)
sim.simulate(sys)

print bb.calcPressureDrop(2.1e6,3.927,0.2503,3,800,0.00003,0.002,10,0.00001,.2,100,50)
print pa.calcPressureDrop(2.1e6,3.927,0.2503,3,800,0.00003,0.002,10,0.00001,.2,100,50)

##h,Al,Ag,Sl,Sg,Si,t_wL,t_wG,t_i=pa.calc_h(800,30,.002,3e-5,5,2,0.0001,.2,0)

##d,Af,Ac,Sl,Si,t_wL,t_i,FE=pa.calc_d(800,30,.002,1e-5,5,10,0.0001,.5,1,100,1e6)
##print d
##print pa.calc_d_min(800,30,.002,1e-5,5.0,5.0,0.0001,.5,-.5,100,1e6)

