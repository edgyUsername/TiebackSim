import sys
from os import path
sys.path.append(path.dirname( path.abspath(__file__)) ) 
import settings


f=settings.System(3000,45,300)
f.addPipe(3000,5,5,5,5)
f.setGeometry()
