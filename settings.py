import math


class Pipe:
	def __init__(self, length, diameter, relative_roughness, thickness, inclination):
		assert type(length)==float or type(length)==int, "The length of the pipe must be a float or integar"
		assert type(diameter) or type(diameter)==int, "The diameter of the pipe must be a float or integar"
		assert type(relative_roughness) or type(relative_roughness)==int, "The relative_roughness of the pipe must be a float or integar"
		assert type(thickness) or type(thickness)==int, "The thickness of the pipe must be a float or integar"
		assert type(inclination) or type(inclination)==int, "The inclination of the pipe must be a float or integar"
		self.l=float(length)
		self.d=float(diameter)
		self.e=float(relative_roughness)
		self.t=float(thickness)
		self.a=float(inclination)
		self.z=float(self.l*math.sin(math.radians(self.a)))
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
	def __init__(self,p_in, massFlow, T_in):
		self.sections=[] #list of nodes, pipes, and legs
		self.simData={'t':[],'x':[],'d':[],'y':[],'i':[],'vars':{'T':[],'P':[],'v':[],'vapQ':[],'r':[]}} #t is time; x is x coord;y is y coord; d is distance; i is index of pipe; T is temp; P iss pressure; v is velocity; vap Q is mass vap quality; r is density
		self.p_in=p_in
		self.T_in=T_in
		self.massFlow=massFlow
		self.__addNode(p_in,T_in,massFlow)
		self.devisions=100
		self.lengths=[0]

	def __addNode(self,p,T, massFlow, prev=None, next=None):
		k=Node()
		k.p=p
		k.massFlow=massFlow
		k.T=T
		self.sections.append(k)

	def addPipe(self,length, diameter, relative_roughness, thickness, inclination):
		#add a pipe to system and connect nodes
		k=self.sections[-1]
		i=len(self.sections)
		p=Pipe(length, diameter, relative_roughness, thickness, inclination)
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
		while d[-1]+delta<tot_len:
			d.append(d[-1]+delta)
		d.append(tot_len)
		x=[0]
		y=[0]
		index=[]
		c=1
		p=self.sections[alloc[0]]
		for j in d[:-1]:
			changePipe=False
			while j>l[c]:
				changePipe=True
				c+=1
			if changePipe:
				p=self.sections[alloc[c-1]]
			deltaX=p.x*delta/p.l
			deltaY=p.z*delta/p.l
			x.append(x[-1]+deltaX)
			y.append(y[-1]+deltaY)
			index.append(alloc[c-1])
		index.append(alloc[-1])
		self.simData['d']=d
		self.simData['x']=x
		self.simData['y']=y
		self.simData['i']=index








		




