import os
import pathlib
from kevtools_common.types import TCID
import json
import warnings
from Bio import SeqIO
import gzip

class Deuterocol1(object):
	def __init__(self, famlist=None, tmdatadir=None, bundle=None, level=3):
		if famlist is None: self.famlist = []
		else: self.famlist = [TCID(fam) for fam in famlist]

		self.backward = False

		self.tmdatadir = pathlib.Path(tmdatadir) if tmdatadir is not None else tmdatadir
		self.bundle = bundle
		self.level = level

		self.pdbclist = []
		self.indices = {}
		self.tcmap = {}
		self.fammap = {}
		self.pdblengths = {}

	def fetch(self):
		pdbclist = []
		tcmap = {}
		fammap = {}
		with open('{}/tcmap.tsv'.format(self.tmdatadir)) as fh:
			for l in fh:
				tcid, strucs = l.strip().split()
				tcid = TCID(tcid)
				for fam in self.famlist: 
					if tcid in fam:
						struclist = strucs.split(',')
						pdbclist.extend(struclist)
						tcmap.update(dict([[pdbc, tcid] for pdbc in struclist]))
						if fam not in fammap: fammap[fam] = []
						fammap[fam].extend(struclist)
						break

		self.tcmap = tcmap
		self.fammap = fammap
		self.pdbclist = pdbclist
		pdbcset = set(pdbclist)


		self.pdblengths = {}

		#legacy tmdata
		if self.tmdatadir.joinpath("opm/ASSIGNMENTS.TSV").exists():
			with open(self.tmdatadir.joinpath("opm/ASSIGNMENTS.TSV")) as fh:
				for l in fh:
					sl = l.strip().split('\t')
					if len(sl) == 0: continue
					elif l.startswith('#'): continue
					elif len(sl) == 1:
						pdbc = sl[0]
						indices = []
					else:
						pdbc, indices = sl
						indices = [[int(x) for x in span.split('-')] for span in indices.split(',')]
					if pdbc in pdbcset:
						self.indices[pdbc] = indices


			for pdbc in self.indices:
				if not self.indices[pdbc]: tmlength = 0
				else: tmlength = sum([span[1] - span[0] + 1 for span in self.indices[pdbc]])

				try: totlength = len(SeqIO.read('{}/sequences/{}.faa'.format(self.tmdatadir, pdbc), 'fasta').seq)
				except FileNotFoundError:
					if self.indices[pdbc]: totlength = self.indices[pdbc][-1][-1] - self.indices[pdbc][0][0] + 1
					else: totlength = None
				self.pdblengths[pdbc] = {'tm':tmlength, 'total':totlength}

		elif self.tmdatadir.joinpath("assignments/opm_tms_extended.json.gz").exists():
			with gzip.open(self.tmdatadir.joinpath("assignments/opm_tms_extended.json.gz"), "rt") as fh:
				self.indices = json.load(fh)

			totlengths = {}
			for record in SeqIO.parse(self.tmdatadir.joinpath("sequences/gapless.faa"), "fasta"):
				totlengths[record.id] = len(record)
			for pdbc in self.indices:
				if not self.indices[pdbc]: tmlength = 0
				else: tmlength = sum([span[1] -span[0] + 1 for span in self.indices[pdbc]])
				self.pdblengths[pdbc] = {"tm":tmlength,
					"total":totlengths.get(pdbc, -1)
				}


	def get_pdblength(self, pdbc, tmrange=None, resirange=None, masktm=0):
		if tmrange is None and resirange is None: 
			return self.pdblengths[pdbc[:4].lower() + pdbc[4:]]

		elif tmrange:
			if mask == 1:
				return sum([self.get_pdblength(pdbc, tms, resirange, masktm=0) for tms in range(tmrange[0], tmrange[1]+1)]) 
			elif mask == -1:
				#this is really bad and inaccurate
				return self.get_pdblength(pdbc, tmrange, resirange, 0) - self.get_pdblength(pdbc, tmrange, resirange, 1)
			else:
				start, stop = tmrange
				indices = self.indices[pdbc][start-1:stop+1-1]
				resirange = [indices[0][0], indices[-1][-1]]
				return self.get_pdblength(pdbc, resirange=resirange)

		elif resirange:
			start, stop = resirange
			found = set()
			with open('{}/pdb/{}.pdb'.format(self.tmdatadir, pdbc.lower()[:4])) as fh:
				for l in fh:
					if l.startswith('ATOM  ') or l.startswith('HETATM'):
						if (pdbc[5] == l[21]) or (pdbc[5] == 'A' and l[21] in 'A '):
							n = int(l[22:26])
							if start <= n <= stop: found.add(n)
			return len(found)

	def get_states(self, pdbc, tmrange=None, resirange=None, indexed=False):
		if tmrange is None and resirange is None:
			return self.get_states(pdbc, resirange=[-1, 9999], indexed=indexed) #as long as we don't have any full structures of titin, this should be ok
		elif tmrange:
			start, stop = tmrange
			indices = self.indices[pdbc][start-1:stop+1-1]
			resirange = [indices[0][0], indices[-1][-1]]
			return self.get_states(pdbc, resirange=[resirange[0], resirange[1]], indexed=indexed)
		elif resirange:
			indices = self.indices[pdbc[:4].lower() + pdbc[4:]]
			start, stop = resirange
			states = {}

			if self.tmdatadir.joinpath('pdb/{}.pdb'.format(pdbc.lower()[:4])).exists(): #legacy
				fh = open(self.tmdatadir.joinpath('pdb/{}.pdb'.format(pdbc.lower()[:4])))
			elif self.tmdatadir.joinpath('structures/pdb/{}.pdb'.format(pdbc.lower()[:4])).exists(): #new
				fh = open(self.tmdatadir.joinpath('structures/pdb/{}.pdb'.format(pdbc.lower()[:4])))
			else: raise FileNotFoundError('Could not find a suitable PDB structure for {}'.format(repr(pdbc)))

			for l in fh:
				if l.startswith('ATOM  ') or l.startswith('HETATM'):
					if (pdbc[5] == l[21]) or (pdbc[5] == 'A' and l[21] in 'A '):
						n = int(l[22:26])
						if start <= n <= stop: 
							states[n] = self._in_spans(n, indices) 

			fh.close()

			if indexed: return states
			else: return [states[i] for i in sorted(states)]
			
	def _in_spans(self, n, spans):
		for start, stop in spans:
			if start <= n <= stop: return True
		return False

	def __contains__(self, thing):
		return self.tcmap.__contains__(thing)


	def __add__(self, other):
		d1sum = Deuterocol1()
		
		d1sum.famlist = self.famlist + other.famlist
		if self.tmdatadir: d1sum.tmdatadir = self.tmdatadir
		else: d1sum.tmdatadir = other.tmdatadir

		if self.bundle: d1sum.bundle = self.bundle
		else: d1sum.bundle = other.bundle

		d1sum.level = max(self.level, other.level)

		d1sum.indices.update(self.indices)
		d1sum.indices.update(other.indices)

		d1sum.tcmap.update(self.tcmap)
		d1sum.tcmap.update(other.tcmap)

		d1sum.fammap.update(self.fammap)
		d1sum.fammap.update(other.fammap)

		d1sum.pdblengths.update(self.pdblengths)
		d1sum.pdblengths.update(other.pdblengths)

		return d1sum

	def write(self, outdir):
		if not os.path.isdir(outdir): os.mkdir(outdir)

		with open('{}/indices.json'.format(outdir), 'w') as fh:
			fh.write(json.dumps(self.indices).replace(']], ', ']],\n'))

		with open('{}/pdblist.json'.format(outdir), 'w') as fh:
			sanitized = dict([[str(k), self.fammap[k]] for k in self.fammap])
			fh.write(json.dumps(sanitized).replace('], ', '],\n'))

		with open('{}/tcmap.json'.format(outdir), 'w') as fh:
			sanitized = dict([[k, str(self.tcmap[k])] for k in self.tcmap])
			fh.write(json.dumps(sanitized).replace(', ', ',\n'))

		with open('{}/pdblengths.json'.format(outdir), 'w') as fh:
			fh.write(json.dumps(self.pdblengths))

		with open('{}/config.json'.format(outdir), 'w') as fh:
			fh.write(json.dumps({
				'famlist':[str(x) for x in self.famlist],
				'tmdatadir':str(self.tmdatadir),
				'bundle':int(self.bundle),
				'level':int(self.level),
			}))

	def read(self, outdir):
		try:
			with open('{}/config.json'.format(outdir)) as fh:
				obj = json.load(fh)
				self.famlist = [TCID(x) for x in obj['famlist']]
				self.backward = False
		except FileNotFoundError: 
			warnings.warn('Could not find config.json. Starting backward-compatible mode')
			self.backward = True

		with open('{}/indices.json'.format(outdir)) as fh:
			self.indices = json.load(fh)

		with open('{}/pdblist.json'.format(outdir)) as fh:
			sanitized = json.load(fh)
			self.fammap = dict([[TCID(k), sanitized[k]] for k in sanitized])

		with open('{}/tcmap.json'.format(outdir)) as fh:
			sanitized = json.load(fh)
			self.tcmap = dict([[k, TCID(sanitized[k])] for k in sanitized])
		
		with open('{}/pdblengths.json'.format(outdir)) as fh:
			self.pdblengths = json.load(fh)

def load_family(famlist, tmdatadir, bundle=None, level=3):
	obj = Deuterocol1(famlist, tmdatadir, bundle=bundle, level=level)
	obj.fetch()
	return obj

def write(d1list, outdir):
	d1sum = Deuterocol1()
	for d1 in d1list: d1sum = d1sum + d1

	d1sum.write(outdir)
