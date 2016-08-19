import sys
import csv
from matplotlib import pyplot as plt
from os import path
sys.path.append(path.dirname( path.abspath(__file__)) ) 


def write_geometry_file(system, disp=False):
	try:
		with open('results/geometry.csv', 'wb') as csvfile:
			fieldnames = ['x (m)', 'y (m)','pipe id']
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
			writer.writeheader()
			for i in system.simData['geometry']:
				writer.writerow({'x (m)':i['x'], 'y (m)':i['y'],'pipe id':i['index']})
	except IOError:
		print "failed to write geometry file. Make sure geometry.csv is not open by other program"

def write_summary_file(system, disp=True):
	if disp:
		x=[system.simData['geometry'][i['geo_index']]['d']/1000. for i in system.simData['vars'] ]
		label={'P':r'$P$ (kPa a)','T':r'$T$ ($^\circ$ C)', 'hold_up':r'$\lambda_{l}$','z':r'$\Delta z$ (m)','x':r'$x$ (km)','liqQ':r'$\beta^{\mathrm{l}}$',
			   'r_l':r'$\rho_{L}$ (kg m$^{-3}$)','r_g':r'$\rho_{G}$ (kg m$^{-3}$)'}
		z=[system.simData['geometry'][i['geo_index']]['y'] for i in system.simData['vars']]
		P=[i['P']/1e3 for i in system.simData['vars']]
		T=[i['T'] for i in system.simData['vars']]
		D=[i['r_l'] for i in system.simData['vars']]
		D2=[i['r_g'] for i in system.simData['vars']]
		hold_up=[i['hold_up'] for i in system.simData['vars']]
		liqQ=[1-i['vapQ'] for i in system.simData['vars']]
		fig=plt.figure()
		plotGeo = fig.add_subplot(321)
		plotP=fig.add_subplot(322)
		plotT = fig.add_subplot(323)
		ploth = fig.add_subplot(324)
		plotV = fig.add_subplot(326)
		plotD= fig.add_subplot(325)

		plotGeo.plot(x,z)
		plotP.plot(x, P)
		plotT.plot(x, T)
		ploth.plot(x, hold_up)
		plotV.plot(x, liqQ)
		plotD.plot(x, D)
		plotD2 = plotD.twinx()
		plotD2.plot(x,D2,'r-')
		plotGeo.set_xlabel(label['x'])
		plotP.set_xlabel(label['x'])
		plotT.set_xlabel(label['x'])
		ploth.set_xlabel(label['x'])
		plotV.set_xlabel(label['x'])
		plotD.set_xlabel(label['x'])
		# plotD2.set_xlabel(label['x'])
		plotGeo.set_ylabel(label['z'])
		plotP.set_ylabel(label['P'])
		plotT.set_ylabel(label['T'])
		ploth.set_ylabel(label['hold_up'])
		plotV.set_ylabel(label['liqQ'])
		plotD.set_ylabel(label['r_l'])
		plotD2.set_ylabel(label['r_g'],color='r')
		for tl in plotD2.get_yticklabels():
			tl.set_color('r')
		plt.tight_layout()
		plt.show()
		print
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