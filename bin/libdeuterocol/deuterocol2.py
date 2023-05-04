import os
import json
import pathlib

class Deuterocol2(object):
	def __init__(self, d1list=None, tmdatadir=None, bundle=None, allow_internal=False, cut_loops=0):
		''' Constructor for Deuterocol2

		bundle (default: None): Size of TMS/loop bundles to use. A value of None bypasses the cutting step, allowing for whole-structure alignments
		cut_loops (default: 0): Loop-cutting behavior. 0 skips cuts, 1 cuts loops away, and -1 cuts things that are not loops away

		'''

		self.d1list = [] if d1list is None else d1list
		self.bundle = bundle
		self.tmdatadir = pathlib.Path(tmdatadir) if tmdatadir is not None else tmdatadir

		if self.tmdatadir is not None:
			if self.tmdatadir.joinpath("pdb").exists():
				self.pdbdir = self.tmdatadir.joinpath("pdb")
			elif self.tmdatadir.joinpath("structures/pdb").exists():
				self.pdbdir = self.tmdatadir.joinpath("structures/pdb")
			else:
				raise FileNotFoundError("Could not find a valid PDB directory in {}".format(self.tmdatadir))

		self.allow_internal = allow_internal
		self.cut_loops = cut_loops

		self.fragments = {}

	def compute(self):
		fragments = {}
		fammaps = {}
		for index, deut in enumerate(self.d1list):
			fragments[index] = []

			fammaps[index] = revfammap(deut.fammap)

			if self.bundle is None:
				for pdbc in deut.indices: 
					if len(deut.indices[pdbc]): label = 'h{}-{}'.format(1, len(deut.indices[pdbc]))
					else: label = 's0'

					fragments[index].append({'pdbc':pdbc, 'span':None, 'label':label})

			elif self.cut_loops == 0:
				#no loop cutting
				for pdbc in deut.indices:
					if len(deut.indices[pdbc]) == 0: continue
					elif len(deut.indices[pdbc]) <= self.bundle: 
						fragments[index].append({'pdbc':pdbc, 'span':[[deut.indices[pdbc][0][0], deut.indices[pdbc][-1][-1]]], 'label':'h{}-{}'.format(1, len(deut.indices[pdbc]))})
					else:
						for i in range(len(deut.indices[pdbc]) - self.bundle + 1):
							fragments[index].append({'pdbc':pdbc, 'span':[[deut.indices[pdbc][i][0], deut.indices[pdbc][i+self.bundle-1][-1]]], 'label':'h{}-{}'.format(i+1, i+self.bundle)})

			elif self.cut_loops == 1:
				#loop cutting
				for pdbc in deut.indices:
					if len(deut.indices[pdbc]) == 0: continue
					elif len(deut.indices[pdbc]) <= self.bundle:
						fragments[index].append({'pdbc':pdbc, 'span':deut.indices[pdbc], 'label':'m{}-{}'.format(1, len(deut.indices[pdbc]))})
					else:
						for i in range(len(deut.indices[pdbc]) - self.bundle + 1):
							fragments[index].append({'pdbc':pdbc, 'span':deut.indices[pdbc][i:i+self.bundle], 'label':'m{}-{}'.format(i+1, i+self.bundle)})

			elif self.cut_loops == -1:
				#loop sparing
				for pdbc in deut.indices:
					if len(deut.indices[pdbc]) == 0: fragments[index].append({'pdbc':pdbc, 'span':None, 'label':'s0'})

					elif len(deut.indices[pdbc]) < (self.bundle + 2):
						spans = []
						spans.append([None, deut.indices[pdbc][0][0]])
						for i in range(len(deut.indices[pdbc])-1):
							spans.append([deut.indices[pdbc][i][-1], deut.indices[pdbc][i+1][0]])
						spans.append([deut.indices[pdbc][-1][0], None])

						fragments[index].append({'pdbc':pdbc, 'span':spans, 'label':'s{}-{}'.format(0, len(deut.indices[pdbc]))})

					else:
						bundlecount = len(deut.indices[pdbc]) - self.bundle + 2
						loops = []
						for i in range(len(deut.indices[pdbc])+1):
							if i == 0: start = None
							else: start = deut.indices[pdbc][i-1][-1]

							if i == (len(deut.indices[pdbc])): stop = None
							else: stop = deut.indices[pdbc][i][0]

							loops.append([start, stop])

						for i in range(bundlecount):
							
							if self.bundle == 1: label = 's{}'.format(i)
							else: label = 's{}-{}'.format(i, i+self.bundle-1)

							fragments[index].append({'pdbc':pdbc, 'span':loops[i:i+self.bundle], 'label':label})
			else: raise ValueError('Unknown loop-cutting mode {}'.format(self.cut_loops))
		#self.fragments = fragments

		self.fragments = {}
		for index in fragments:
			self.fragments[index] = {}
			for fragment in fragments[index]:
				pdbc = "{}_{}".format(fragment['pdbc'][:4].upper(), fragment['pdbc'][-1])
				try: fam = fammaps[index][pdbc]
				except KeyError: continue
				if fam not in self.fragments[index]: self.fragments[index][fam] = []
				self.fragments[index][fam].append(fragment)

	def get_alignments(self):
		for index1 in self.fragments:
			for index2 in self.fragments:
				if index1 >= index2: continue
				for fam1 in self.fragments[index1]:
					for fam2 in self.fragments[index2]:
						for frag1 in self.fragments[index1][fam1]:
							for frag2 in self.fragments[index2][fam2]:
								yield frag1, frag2

	def get_fragments(self):
		for index in self.fragments:
			for fam in self.fragments[index]:
				for frag in self.fragments[index][fam]: yield frag

	def write(self, outdir):
		if not os.path.isdir(outdir): os.mkdir(outdir)

		for index1 in self.fragments:
			for index2 in self.fragments:
				if index1 >= index2: continue
				for fam1 in self.fragments[index1]:
					for fam2 in self.fragments[index2]:
						pairpath = '{}/{}_vs_{}'.format(outdir, fam1, fam2)
						if not os.path.isdir(pairpath): os.mkdir(pairpath)

						if not os.path.isdir('{}/config'.format(pairpath)): os.mkdir('{}/config'.format(pairpath))

						with open('{}/config/align_me.json'.format(pairpath), 'w') as fh:
							obj = {
								str(fam1): [frag['pdbc'] for frag in self.fragments[index1][fam1]],
								str(fam2): [frag['pdbc'] for frag in self.fragments[index2][fam2]],
							}
							fh.write(json.dumps(obj))

						with open('{}/config/indices.json'.format(pairpath), 'w') as fh:
						#FIXME: Why didn't I make this an actual JSON the first time around exactly???
							[fh.write('{pdbc}\t{span}\n'.format(**fragment)) for fragment in self.fragments[index1][fam1]]
							[fh.write('{pdbc}\t{span}\n'.format(**fragment)) for fragment in self.fragments[index2][fam2]]

						with open('{}/config/config.json'.format(pairpath), 'w') as fh:
							obj = {
								'allow_internal':self.allow_internal,
								'bundle':self.bundle,
								'cut_loops':self.cut_loops,
							}
							fh.write(json.dumps(obj))

						with open('{}/config/agenda.json'.format(pairpath), 'w') as fh: 
							#FIXME: Also no true JSON! Whyyyyyyy?
							for frag1 in self.fragments[index1][fam1]:
								for frag2 in self.fragments[index2][fam2]:
									obj = {}
									obj['query'] = frag1['pdbc'][:4]
									obj['subject'] = frag2['pdbc'][:4]

									obj['qchain'] = frag1['pdbc'][-1]
									obj['schain'] = frag2['pdbc'][-1]

									obj['qlabel'] = frag1['label']
									obj['slabel'] = frag2['label']

									if frag1['span'] is None: obj['qspan'] = [1, self.d1list[index1].pdblengths[frag1['pdbc']]]
									else: obj['qspan'] = [frag1['span'][0][0], frag1['span'][-1][-1]]
									if frag2['span'] is None: obj['sspan'] = [1, self.d1list[index2].pdblengths[frag2['pdbc']]]
									else: obj['sspan'] = [frag2['span'][0][0], frag2['span'][-1][-1]]
									
									obj['qspan'] = frag1['span']
									obj['sspan'] = frag2['span']
									obj['bundle'] = self.bundle
									obj['qindices'] = frag1['span']
									obj['sindices'] = frag2['span']
									obj['name'] = '{}_{}_vs_{}_{}'.format(frag1['pdbc'], frag1['label'], frag2['pdbc'], frag2['label'])

									fh.write(json.dumps(obj) + '\n')
							
						with open('{}/config/tcmap.json'.format(pairpath), 'w') as fh: 
							obj = {}
							for deut in self.d1list:
								for pdbc in deut.tcmap:
									obj[pdbc] = str(deut.tcmap[pdbc])
							fh.write(json.dumps(obj))


		#for aln in self.get_alignments(): print(aln)

		#if not os.path.isdir('{}/config'.format(outdir)): os.mkdir('{}/config'.format(outdir))
	def get_pdblength(self, pdbc, tmrange=None, resirange=None, masktm=None):
		for library in self.d1list:
			if pdbc in library: 
				return library.get_pdblength(pdbc, tmrange=tmrange, resirange=resirange, masktm=masktm)
			elif pdbc[:4].upper() + pdbc[4:] in library: 
				return library.get_pdblength(pdbc[:4].upper() + pdbc[4:], tmrange=tmrange, resirange=resirange, masktm=masktm)
		raise KeyError("Could not find {} in library".format(pdbc))

	def get_states(self, pdbc, tmrange=None, resirange=None, indexed=False, masktm=None):
		for library in self.d1list:
			if pdbc in library:
				return library.get_states(pdbc, tmrange=tmrange, resirange=resirange, indexed=indexed)
			elif pdbc[:4].upper() + pdbc[4:] in library: 
				return library.get_states(pdbc[:4].upper() + pdbc[4:], tmrange=tmrange, resirange=resirange, indexed=indexed)

				#return library.get_pdblength(pdbc, tmrange=tmrange, resirange=resirange, masktm=masktm)
				#return library.get_pdblength(pdbc[:4].upper() + pdbc[4:], tmrange=tmrange, resirange=resirange, masktm=masktm)
		raise KeyError("Could not find {} in library".format(pdbc))

def revfammap(fammap):
	revmap = {}
	for fam in fammap:
		for pdbc in fammap[fam]:
			revmap[pdbc] = fam
	return revmap
