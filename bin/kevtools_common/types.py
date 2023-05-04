import re
import json

class TCID:
	def __init__(self, somestring=None, somelist=None):
		self.iterable = [None, None, None, None, None, None]

		if somestring: self.parse_string(somestring)
		elif somelist:
			for i, x in enumerate(somelist): self.iterable[i] = x
			self.validate()

	def parse_string(self, somestring):
		vals = re.split('[\.\-]', somestring.strip())
		for i, x in enumerate(vals[:6]): self.iterable[i] = x
		self.validate()

	def validate(self):
		for i, x in enumerate(self.iterable):
			if i == 0 and x is not None: assert re.match('^[0-9]$', x), 'Invalid class id "{}" for {} (1st field)'.format(x, self)
			elif i == 1 and x is not None: assert re.match('^[A-Z]$', x), 'Invalid subclass id "{}" for {} (2nd field)'.format(x, self)
			elif i == 2 and x is not None: assert re.match('^[0-9]+$', x), 'Invalid family id "{}" for {} (3rd field)'.format(x, self)
			elif i == 3 and x is not None: assert re.match('^[0-9]+$', x), 'Invalid subfamily id "{}" for {} (4th field)'.format(x, self)
			elif i == 4 and x is not None: assert re.match('^[0-9]+$', x), 'Invalid system id "{}" for {} (5th field)'.format(x, self)
			elif i == 5 and x is not None: assert True, 'Invalid external accession "{}" for {} (6th field)'.format(x, self)

		must_be_defined = False
		for i, x in enumerate(self.iterable[::-1]):
			if must_be_defined: assert x is not None, 'Invalid TCID: Undefined parent field'
			elif x is not None: must_be_defined = True

	@staticmethod
	def from_string(somestring):
		tcobj = TCID(somestring=somestring)
		return tcobj

	@staticmethod
	def from_list(somelist):
		tcobj = TCID(somelist=somelist)
		return tcobj

	def __iter__(self): return iter(self.iterable[:6])

	def __len__(self): return sum([x is not None for x in self])

	def __contains__(self, other):
		for x1, x2 in zip(self, other):
			if x1 is None and x2 is None: pass
			elif x1 is None and x2 is not None: pass
			elif x1 is not None and x2 is not None:
				if x1 != x2: return False
			elif x1 is not None and x2 is None: return False
		return True

	def __getitem__(self, index): 
		if type(index) is slice: return TCID.from_list(self.iterable[index])
		elif type(index) is int: 
			#if index in (0, 2, 3, 4): 
			#	if self.iterable[index] is None: return None
			#	else: return int(self.iterable[index])
			#else: 
			return self.iterable[index]
		else: raise TypeError

	def __repr__(self):
		return 'TCID("{}")'.format(self.to_string())

	def to_string(self):
		s = ''
		for i, x in enumerate(self):
			if x is None: break

			if i == 0: pass
			elif 0 < i < 5: s += '.'
			elif i == 5: s += '-'
			else: s += '.'

			s += x
		return s

	def __str__(self):
		return self.to_string()

	def __eq__(self, other):
		for x1, x2 in zip(self, other):
			if x1 != x2: return False
		return True

	def __ne__(self, other):
		return not self.__eq__(other)

	def __lt__(self, other):
		if self.__eq__(other): return False

		if not isinstance(other, TCID): return False

		for i in range(6):
			if self[i] is None and other[i] is None: continue
			elif (self[i] is None) ^ (other[i] is None): return str(self) < str(other)

			if i == 0:
				if int(self[0]) > int(other[0]): return False
				elif int(self[0]) < int(other[0]): return True
			elif i == 1:
				if self[1] > other[1]: return False
				elif self[1] < other[1]: return True
			elif i == 2:
				if int(self[2]) > int(other[2]): return False
				elif int(self[2]) < int(other[2]): return True
			elif i == 3:
				if int(self[3]) > int(other[3]): return False
				elif int(self[3]) < int(other[3]): return True
			elif i == 4:
				if int(self[4]) > int(other[4]): return False
				elif int(self[4]) < int(other[4]): return True
			elif i == 5:
				if self[5] > other[5]: return False
				elif self[5] < other[5]: return True

		return True

	def __gt__(self, other):
		if self.__eq__(other): return False
		return not self.__lt__(other)

	def __hash__(self):
		return hash(str(self)) ^ hash(tuple(self.iterable))

class TCIDCollection(object):
	def __init__(self, iterable=None):
		self.iterable = set()
	
		if iterable:
			for elem in iterable:
				if type(elem) is str: self.iterable.add(TCID(somestring=elem))
				elif type(elem) is TCID: self.iterable.add(elem)
				else: 
					print(type(elem))
					raise TypeError('iterable must contain TCIDs')


	def __len__(self): return len(self.iterable)
	def add_tcid(self, tcid):
		if type(tcid) is str: tcid = TCID(somestring=tcid)
		self.iterable.add(tcid)
	def remove_tcid(self, tcid):
		if type(tcid) is str: tcid = TCID(somestring=tcid)
		self.iterable.discard(tcid)

	def __iter__(self): return iter(self.iterable)

	def __repr__(self):
		out = '{'

		for tcid in self: out += repr(tcid) + ', '
		if out.endswith(', '): out = out[:-2]

		out += '}'
		return out

	def search(self, query):
		if type(query) is str: query = TCID(somestring=query)
		matches = TCIDCollection()
		for tcid in self:
			if tcid in query: matches.add_tcid(tcid)
		return matches

class GODB(object):
	def __init__(self):
		self.go2acc = {}
		self.acc2go = {}
		self.golabels = {}

	def serialize(self):
		return json.dumps({
			'go2acc':self.go2acc,
			'acc2go':self.acc2go,
			'golabels':self.golabels,
		}, indent=1)
	@staticmethod
	def deserialize(obj):
		newdb = GODB()
		newdb.go2acc = obj['go2acc']
		newdb.golabels= obj['golabels']
		newdb.acc2go = obj['acc2go']
		return newdb

	def add(self, go=None, acc=None, label=None):
		if go not in self.go2acc: self.go2acc[go] = []
		if acc not in self.acc2go: self.acc2go[acc] = []
		self.go2acc[go].append(acc)
		self.acc2go[acc].append(go)
		self.golabels[go] = label
