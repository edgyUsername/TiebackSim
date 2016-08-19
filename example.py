import settings
import data_writer
from simulation import simulator as sim

sysParams= {
    'P0':24.1e5,   #Pa abs
    'T0':50,        #deg C
    # 'Ql':1000.,      #sm3/day
    'm':14,              #kg/s
    'T_amb':4,
    'water_cut':None,
    'sysVars':[
        {
        'type':'Pipe',  #Pipe or Leg
        'params':(
            10000,       #length (m)
            .343,      #diameter (m)
            0.0000254,  #roughness (m)
            0.0127,     #thickness (m)
            -0.01,        #inclination (deg)
            10        #U (W/m2/K)
            )
        },{
        'type':'Pipe',  #Pipe or Leg
        'params':(
            200,       #length (m)
            0.343,      #diameter (m)
            0.0000254,  #roughness (m)
            0.0127,     #thickness (m)
            90.,        #inclination (deg)
            10        #U (W/m2/K)
            )
        }
    ],
    'comp':[
        (36.5,'methane'),
        (4.4,'ethane'),
        (2.6,'propane'),
        (.63,'ibutane'),
        (.13, 'butane'),
        (.67, 'ipentane'),
        (.83, 'pentane'),
        (2.7, 'hexane'),
        (51.54,'octane')
        # (78.03,'methane'),
        # (9.7,'ethane'),
        # (4.82,'propane'),
        # (.6,'isobutane'),
        # (1.33,'butane'),
        # (.3,'isopentane'),
        # (.3,'pentane'),
        # (3.22,'octane'),
        # (.14,'he'),
        # (.013,'hydrogen'),
        # (4.55,'nitrogen'),
        # (78.03,'carbondioxide'),
    ]
}
"""############### Unit Conversions #####################
1bpd=1.840130787037037e-06		                	m3/s 
1scfd=3.2774128000000003e-07		                m3/s
1 lb/ft3=16.018					                  	kg/m3
1"=0.0254			    		                  	m
1psi=6894.757						                Pa
######################################################"""
sys=settings.newSystem(sysParams)

sys=sim.forward_steady_itterate(sys,VLE='pr',flow='p_a',tol_P=1,devs=1000)

data_writer.write_summary_file(sys)

