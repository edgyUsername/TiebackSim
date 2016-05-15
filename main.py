import sys
from os import path
sys.path.append(path.dirname( path.abspath(__file__)) ) 
import settings




sysParams= {
    'P0':2410000,   #Pa abs
    'T0':50,        #deg C
    'Ql':3280,      #sm3/day
    'sysVars':[
        {
        'type':'Pipe',  #Pipe or Leg
        'params':(
            5000,       #length (m) 
            0.343,      #diameter (m)
            0.0000254,  #roughness (m)
            0.0127,     #thickness (m)
            0.0         #inclination (deg)
            )
        },{
        'type':'Pipe',  
        'params':(
            5000,       
            0.343,      
            0.0000254,  
            0.0127, 
            -0.1718875964         
            )
        },{
        'type':'Pipe',  
        'params':(
            200,       
            0.343,      
            0.0000254,  
            0.0127, 
            90         
            )
        }
    ]
}


sys=settings.newSystem(sysParams)
sys.printGeometry()
