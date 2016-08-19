import math
import simulation.thermodynamics.PVT as PVT
from simulation.thermodynamics import setup_comp_data as configPVT
ln=math.log
log=math.log10
exp=math.exp
atan=math.atan
sqrt=math.sqrt
class Pipe:
	def __init__(self, length, diameter, roughness, thickness, inclination, U):
		assert type(length)==float or type(length)==int, "The length of the pipe must be a float or integar"
		assert type(diameter)==float or type(diameter)==int, "The diameter of the pipe must be a float or integar"
		assert type(roughness)==float or type(roughness)==int, "The roughness of the pipe must be a float or integar"
		assert type(thickness)==float or type(thickness)==int, "The thickness of the pipe must be a float or integar"
		assert type(inclination)==float or type(inclination)==int, "The inclination of the pipe must be a float or integar"
		assert type(U)==float or type(U)==int, "The overall heat trasfer coefficient, U, must be a float or integar"
		self.l=float(length)
		self.d=float(diameter)
		self.e=float(roughness)
		self.t=float(thickness)
		self.a=float(inclination)
		self.z=float(self.l*math.sin(math.radians(self.a)))
		self.U=float(U)
		self.cosA=math.cos(math.radians(self.a))
		self.sinA=math.sin(math.radians(self.a))
		self.x=length*self.cosA
		self.pIn=None
		self.pOut=None
		self.massFlow=None
		self.TIn=None
		self.TOut=None

class Leg:
	def __init__(self, pDrop, absolute):
		assert type(pDrop)==float or type(pDrop)==int, "The pressure drop of the leg must be a float or integar"
		assert type(absolute)==bool, "absolute must be a boolean. True for absolute, False for relative"
		self.pIn=None
		self.pOut=None
		self.pDrop=pDrop
		self.abs=absolute
		self.massFlow=None
		self.TIn=None
		self.TOut=None

class Node:
	def __init__(self, prev=None, next=None):
		if prev:
			assert type(prev)==int, "prev must be int denoting index of previous pipe or leg"
		if next:
			assert type(next)==int, "next must be int denoting index of next pipe or leg"
		self.prev=prev
		self.next=next
		self.p=None
		self.massFlow=None
		self.T=None

class System:
	def __init__(self,p_in,T_in,T_amb, PVT, pvtType,p_out=None, massFlow=None):
		self.sections=[] #list of nodes, pipes, and legs
		self.simData={
			'vars':[
				
				{
					'geo_index':None,'T':None,'P':None,'v_sg':None,'v_sl':None,'vapQ':None,'r_g':None,'r_l':None,'m_g':None,'m_l':None,'pattern':None,'hold_up':None,'Cp_l':None,'Cp_g':None
				}
				
			],'geometry':[
				{'index':None,'x':None,'y':None,'d':None}
			]
		}#t is time; x is x coord;y is y coord; d is distance; i is index of pipe; T is temp; P iss pressure; v is velocity; vap Q is mass vap quality; r is density
		self.transSimData=[{'t':0,'simData':{}}]
		self.pvtType=pvtType 
		self.pvt=PVT
		self.p_in=p_in
		self.T_in=T_in
		self.p_out=p_out
		self.massFlow=massFlow
		self.__addNode(p_in,T_in,massFlow)
		self.devisions=float(100000)
		self.lengths=[0]
		self.T_amb=T_amb if T_amb else T_in

	def __addNode(self,p,T, massFlow, prev=None, next=None):
		k=Node()
		k.p=p
		k.massFlow=massFlow
		k.T=T
		self.sections.append(k)

	def addPipe(self,length, diameter, roughness, thickness, inclination, U):
		#add a pipe to system and connect nodes
		k=self.sections[-1]
		i=len(self.sections)
		p=Pipe(length, diameter, roughness, thickness, inclination, U)
		p.massFlow=k.massFlow
		p.pIn=k.p
		p.TIn=k.T
		k.next=i
		j=Node(prev=i)
		j.p=p.pOut
		j.T=p.TOut
		j.massFlow=p.massFlow
		self.sections.append(p)
		self.sections.append(j)

		self.lengths.append(self.lengths[-1]+length)

	def addLeg(self,pDrop, absolute):
		#add a leg to the system as well as connect nodes
		k=self.sections[-1]
		i=len(self.sections)
		p=Leg(pDrop, absolute)
		p.massFlow=k.massFlow
		p.pIn=k.p
		p.TIn=k.T
		k.next=i
		j=Node(prev=i)
		j.p=p.pOut
		j.T=p.TOut
		j.massFlow=p.massFlow
		self.sections.append(p)
		self.sections.append(j)

	def setGeometry(self):
		#set x,y coordinates and corresponding pipe index
		tot_len=self.lengths[-1]
		delta=tot_len/self.devisions
		l=self.lengths
		c=0
		alloc=[]
		for i in self.sections:
			if isinstance(i,Pipe):
				alloc.append(c)
			c+=1
		d=[0]
		c=1
		while d[-1]+delta<tot_len:
			changePipe=False
			while d[-1]+delta>l[c]:
				c+=1
				changePipe=True
			if changePipe and len(l)>c-1:
				d.append(l[c-1])
			d.append(d[-1]+delta)
		d.append(tot_len)
		x=[0]
		y=[0]
		index=[]
		c=1
		lc=0
		p=self.sections[alloc[0]]
		for j in d[:-1]:
			changePipe=False
			while d[lc+1]>l[c]:
				changePipe=True
				c+=1
			if changePipe:
				p=self.sections[alloc[c-1]]
			deltaX=p.x*(d[lc+1]-j)/p.l
			deltaY=p.z*(d[lc+1]-j)/p.l
			x.append(x[-1]+deltaX)
			y.append(y[-1]+deltaY)
			index.append(alloc[c-1])
			lc+=1
		index.append(alloc[-1])
		c=0
		geometry=[]
		for j in d:
			geometry.append({'d':j,'index':index[c],'x':x[c],'y':y[c]})
			c+=1
		self.simData['geometry']=geometry
	def printGeometry(self):
		print 'x\ty'
		c=0
		for i in self.simData['geometry']:
			print i['x'],'\t',i['y']
			c+=1


def newSystem(systemData):
	keys=systemData.keys()
	assert 'P0' in keys, "'P0' variable is missing"
	assert 'T0' in keys, "'T0' variable is missing"
	assert 'Ql' in keys or 'm' in keys, "'Ql' or 'm' variable is missing"
	assert 'sysVars' in keys, "'sysVars' variable is missing"
	assert 'PVT' in keys or 'comp' in keys,"'PVT'(black oil) or 'comp'(compositional) variables missing"
	pvtType=None
	if 'comp' in keys:
		pvtType='comp'
		pvt=configPVT.get_props_and_equalize(systemData['comp'],wc=systemData['water_cut'])

		fractions=systemData.get('fractions',{})
		comp_tot=sum([i[0] for i in systemData['comp']])
		for f in fractions:
			x=fractions[f]['comp']/comp_tot
			name=f
			v_c=0.3456+2.4224e-4*fractions[f]['mm']-.4443*fractions[f]['sg']+1.131e-3*fractions[f]['mm']*fractions[f]['sg']
			t_c=338+202*log(fractions[f]['mm']-71.2)+(1361*log(fractions[f]['mm'])-2111)*log(fractions[f]['sg'])
			p_c=81.91-29.7*log(fractions[f]['mm']-61.1)+(fractions[f]['sg']-.8)*(159.9-58.7*log(fractions[f]['mm']-53.7))
			P_vap=exp(ln(p_c/1.01325)*(fractions[f]['t_b']/(t_c-fractions[f]['t_b']))*(1-1/.7))
			hf=-0.15600003121561404*fractions[f]['mm']-208732.21599644143
			k={j:0 for j in pvt}
			for j in pvt:
				pvt[j]['k'][name]=0
				k[name]=0
			Kw=((1.8*fractions[f]['t_b'])**(1/3.))/fractions[f]['sg']
			a=1.4651+.2302*Kw
			b=.306469-.16734*fractions[f]['sg']
			c=.001467-.000551*fractions[f]['sg']
			acc=-log(P_vap)-1.
			C_l=lambda C,T,t_c: a*(b+c*(T+273.15))
			h_l=lambda C,T,t_c: a*(b*(T+273.15)+c*((T+273.15)**2)/2)-a*(b*(298)+c*((298)**2)/2)
			C_ig=lambda C,T: a*(b+c*(T+273.15))-8.314*(1.586+.49/(1-(T+273.15)/t_c)+acc*(4.2775+(6.3/((1-(T+273.15)/t_c)*(T+273.15)/t_c))*abs(1-(T+273.15)/t_c)**(1+1/3.)+.4355/(1-(T+273.15)/t_c)))
			h_ig= lambda C,T: (a*(b+c*(T+273.15))-8.314*(1.586+.49/(1-(T+273.15)/t_c)+acc*(4.2775+(6.3/((1-(T+273.15)/t_c)*(T+273.15)/t_c))*abs(1-(T+273.15)/t_c)**(1+1/3.)+.4355/(1-(T+273.15)/t_c)))-
							   a*(b+c*(298))-8.314*(1.586+.49/(1-(298)/t_c)+acc*(4.2775+(6.3/((1-(298)/t_c)*(298)/t_c))*abs(1-(298)/t_c)**(1+1/3.)+.4355/(1-(298)/t_c))))
			pvt[name]={'comp':x,'k':k,'params':{'_id':name,'heat_l':{'h':h_l,'C':C_l},'heat_g':{'h':h_ig,'C':C_ig},'v_c':v_c,
					   'p_c':p_c,'t_c':t_c,'cl':[],'cg':[],'acc':acc,'hf':hf,'mm':fractions[f]['mm']/1e3}}
		comp_tot=sum([pvt[i]['comp'] for i in pvt])
		for i in pvt:
			pvt[i]['comp']=pvt[i]['comp']/comp_tot
		if 'm' in keys:
			m=systemData['m']
		else:
			data_l=PVT.get_props(pvt,systemData['T0'],systemData['P0'])
			density=data_l['r_l']
			liqQ = 1-data_l['vapQ']
			m=systemData['Ql']*density/(liqQ*3600*24)
		
	elif 'PVT' in keys:
		pvt=systemData['PVT']
		pvtType='bo'
		if 'm' in keys:
			m=systemData['m']
		else:
			m=(systemData['Ql']*pvt['r_l']/(1-pvt['vapQ']))/86400		#86400 sec/day
	else:
		m=None
	sys=System(systemData['P0'],systemData['T0'],systemData['T_amb'], pvt,pvtType,massFlow=m)
	for elem in systemData['sysVars']:
		if elem['type']=='Pipe':
			sys.addPipe(*elem['params'])
		elif elem['type']=='Leg':
			sys.addLeg(*elem['params'])
		else:
			assert False, "sysVars element type must be 'Pipe' of 'Leg'"

	sys.setGeometry()
	return sys


def extrapolate():
	"""
	returns list of all values between two indices
	"""	




