import simulation.simulator as sim
import settings
import math
"""
1bpd=1.840130787037037e-06		                	m3/s
1scfd=3.2774128000000003e-07		                m3/s
1 lb/ft3=16.018					                  	kg/m3
1"=0.0254			    		                  	m
1psi=6894.757						                Pa
"""
area=math.pi*.25*(2.047*0.0254)**2
patterns=[]
flows=[]
P=[]
ul=[]
ug=[]
for liq in range(1,100):
    l=10**(liq*5./100.-2)*area*62.4*16.018
    for gas in range(1,100):
        g=10**(gas*5./100.-2)*area*.08*16.018
        w=(l/.018)/((l/.018)+(g/(.79*.028+.21*.032)))
        sysParams = {
            'P0': 1.013e5,  # Pa abs
            'T0': 25,  # deg C
            'm': l + g,  # kg/s
            'T_amb': 25,
            'sysVars': [
                {
                    'type': 'Pipe',  # Pipe or Leg
                    'params': (
                        .03,  # length (m)
                        2.047 * 0.0254,  # diameter (m)
                        .00015 * 0.0254 * 12,  # roughness (m)
                        0.0127,  # thickness (m)
                        0.,  # inclination (deg)
                        10  # U (W/m2/K)
                    )
                }],
            'comp': [
                (w, 'water'),
                (.21 * (1 - w), 'oxygen'),
                (.79 * (1 - w), 'nitrogen'),
            ]
        }
        sys = settings.newSystem(sysParams)
        try:
            sys= sim.forward_steady_itterate(sys, VLE='pr', flow='p_a',devs=1)
            if not sys.simData['vars'][-2]['pattern'] in patterns:
                patterns.append(sys.simData['vars'][-2]['pattern'])
            flows.append(patterns.index(sys.simData['vars'][-2]['pattern']))
            P.append(sys.simData['vars'][-1]['P'])
            ul.append(10**(liq*5./100.-2))
            ug.append(10 ** (gas * 5. / 100. - 2))
        except:
            pass
pass