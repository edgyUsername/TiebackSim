import csv
import sys
from os import path
sys.path.append( path.dirname( path.dirname(path.dirname( path.abspath(__file__)) ) ) )
#import settings
def all_crit_properties():
	"""returns dict of all properties from Perrys_crit.csv"""
	properties={}
	with open(sys.path[-1]+'/databases/Perrys_crit.csv','rb') as csvfile:
		prop_reader=csv.DictReader(csvfile)
		for row in prop_reader:
			if row['mm'] =="":
				pass
			else:
				names=[row['name'].replace(' ','').replace(',','').replace('-','').lower(),]
				for i in row['alt_name'].split(','):
					names.append(i.replace(' ','').replace('-','').lower()) if not i.replace(' ','').replace('-','') =="" else None
				for n in names:
					properties[n]={
					'p_c':float(row['p_c']),				#bara
					't_c':float(row['t_c']),				#K
					'v_c':float(row['v_c']),				#m3/kmol
					'acc':float(row['acc']),				#accentric factor
					'mm':float(row['mm']),					#kg/mol
					'_id':row['name'].replace(' ','').replace(',','').replace('-','').lower()
					}
	return properties

def get_props_and_equalize(compData):
	"""compData is list with tuples (feed comp, comp name)"""
	ref=all_crit_properties()
	total_moles=0
	for i in compData:
		total_moles+=float(i[0])
	pvtParams={}
	for i in compData:
		try:
			compParams=ref[i[1].replace(',','').replace('-','').replace(' ','').lower()]
		except:
			print '"%s" not available in library or missing some params. You can ammend the Perrys_crit.csv file in the "databases" folder' %i[1]
			raise
		if compParams['_id'] in pvtParams:
			assert False, "components must be unique!"
		pvtParams[compParams['_id']]={'params':compParams,'comp':float(i[0])/total_moles}
	return pvtParams
		
