import sys
import csv
from os import path
sys.path.append(path.dirname( path.abspath(__file__)) ) 
import settings

def write_geometry_file(system):
	try:
		with open('results/geometry.csv', 'wb') as csvfile:
			fieldnames = ['x (m)', 'y (m)','pipe id']
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
			writer.writeheader()
			for i in system.simData['geometry']:
				writer.writerow({'x (m)':i['x'], 'y (m)':i['y'],'pipe id':i['index']})
	except IOError:
		print "failed to write geometry file. Make sure geometry.csv is not open by other program"

def write_summary_file(system):
	try:
		with open('results/summary.csv', 'wb') as csvfile:
			fieldnames = ['d (m)', 'T (deg C)', 'P (Pa)','vg (m/s)','vl (m/s)','vapQ','density g (kg/m3)','density l (kg/m3)','visc g (Pa s)','visc l (Pa s)','Cp g (J/kg/K)','Cp l (J/kg/K)','liquid hold up','flow pattern']
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
			writer.writeheader()
			for i in system.simData['vars']:
				writer.writerow({'d (m)':system.simData['geometry'][i['geo_index']]['d'], 'T (deg C)':i['T'], 'P (Pa)':i['P'],'vg (m/s)':i['v_g'],'vl (m/s)':i['v_l'],'vapQ':i['vapQ'],'density g (kg/m3)':i['r_g'],'density l (kg/m3)':i['r_l'],'visc g (Pa s)':i['m_g'],'visc l (Pa s)':i['m_l'],'Cp g (J/kg/K)':i['Cp_g'],'Cp l (J/kg/K)':i['Cp_l'],'liquid hold up':i['hold_up'],'flow pattern':i['pattern']})
	except IOError:
		print "failed to write summary file. Make sure summary.csv is not open by other program"
def write_pressure_profile(system):
	try:
		with open('results/pressure profile.csv', 'wb') as csvfile:
			fieldnames = ['d (m)','P (Pa)']
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
			writer.writeheader()
			for i in system.simData['vars']:
				writer.writerow({'d (m)':system.simData['geometry'][i['geo_index']]['d'], 'P (Pa)':i['P']})
	except IOError:
		print "failed to write pressure profile file. Make sure pressure profile.csv is not open by other program"

def write_VLE(data):
	try:
		with open(sys.path[-1]+'/results/VLE.csv','wb') as csvfile:
			fieldnames=['T (K)','P (kPa)','Vapor fract']
			writer= csv.DictWriter(csvfile,fieldnames=fieldnames)
			writer.writeheader()
			for i in data:
				writer.writerow({'P (kPa)':i['P']/float(1000),'T (K)':i['T']+273,'Vapor fract':i['V']})
	except IOError:
		print "failed to write VLE file. Make sure VLE.csv is not open by other program"

def write_VLE2(data):
	try:
		with open(sys.path[-1]+'/results/VLEg.csv','wb') as csvfile:
			fieldnames=['T (K)','P (kPa)','Vapor fract']
			writer= csv.DictWriter(csvfile,fieldnames=fieldnames)
			writer.writeheader()
			for i in data:
				writer.writerow({'P (kPa)':i['P']/float(1000),'T (K)':i['T']+273,'Vapor fract':i['V']})
	except IOError:
		print "failed to write VLE file. Make sure VLE.csv is not open by other program"

def write_VLE3(data):
	try:
		with open(sys.path[-1]+'/results/VLEl.csv','wb') as csvfile:
			fieldnames=['T (K)','P (kPa)','Vapor fract']
			writer= csv.DictWriter(csvfile,fieldnames=fieldnames)
			writer.writeheader()
			for i in data:
				writer.writerow({'P (kPa)':i['P']/float(1000),'T (K)':i['T']+273,'Vapor fract':i['V']})
	except IOError:
		print "failed to write VLE file. Make sure VLE.csv is not open by other program"