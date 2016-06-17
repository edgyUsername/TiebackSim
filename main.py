import sys
from os import path
sys.path.append(path.dirname( path.abspath(__file__)) ) 
import settings
import data_writer
from simulation.fluid_flow import b_and_b as bb
from simulation.fluid_flow import p_and_a as pa
from simulation import simulator as sim
from simulation.thermodynamics import setup_comp_data as PVT

sysParams= {
    'P0':2e5,   #Pa abs
    'T0':25,        #deg C
    #'Ql':317.9746,      #sm3/day
    'm':2,              #kg/s
    'T_amb':4,
    'sysVars':[
        {
        'type':'Pipe',  #Pipe or Leg
        'params':(
            10000,       #length (m) 
            0.124,      #diameter (m)
            0.0000254,  #roughness (m)
            0.012,     #thickness (m)
            0.0,        #inclination (deg)
            10.0        #U (W/m2/K)
            )
        }
    ],
    'comp':[
        (21,'oxygen'),
        (79,'nitrogen'),
    ],
    'PVT':{								#black oil pvt data at wellhead
        'r_g':1,						#density of gas (kg/m3)
        'r_l':998.0841,						#density of liquid
        'm_g':0.0000131,					#viscosity of gas (Pa s)
        'm_l':0.00109391,						#viscosity of liquid
        'vapQ':	0.0,			                        #mass vapour quality 
        'Cp_l':77.76696/0.018,
        'Cp_g':0.0
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
#data_writer.write_geometry_file(sys)
#sys.printGeometry()
#sim.simulate(sys,mode='devision',steps=300,devMultiplier=5)

sys,itterations=sim.forward_steady_calculate(sys,VLE='ideal',flow='p_a',tol_P=1)
##sys,itterations=sim.forward_steady_calculate(sys,flow='p_a')
print itterations
data_writer.write_summary_file(sys)
data_writer.write_pressure_profile(sys)
##print bb.calcPressureDrop(2.1e6,3.927,0.2503,3,800,0.00003,0.002,10,0.00001,.2,100,50)
##print pa.calcPressureDrop(2.1e6,3.927,0.2503,3,800,0.00003,0.002,10,0.00001,.2,100,50)

##h,Al,Ag,Sl,Sg,Si,t_wL,t_wG,t_i=pa.calc_h(800,30,.002,3e-5,5,2,0.0001,.2,0)

##d,Af,Ac,Sl,Si,t_wL,t_i,FE=pa.calc_d(800,30,.002,1e-5,5,10,0.0001,.5,1,100,1e6)
##print d
##print pa.calc_d_min(800,30,.002,1e-5,5.0,5.0,0.0001,.5,-.5,100,1e6)

