#!/usr/bin/env python2

from __future__ import print_function, division, generators

import argparse, json, os, subprocess, re, sys, time


VERBOSITY = 1
def info(*things):
	print('[INFO]:', *things, file=sys.stderr)
def warn(*things):
	print('[WARNING]:', *things, file=sys.stderr)

class Alignment(object):
	def __init__(self):
		self.query = ''
		self.subject = ''
		self.queryfn = ''
		self.subjectfn = ''
		self.qhel = []
		self.shel = []
		self.qaligned = []
		self.saligned = []
		#TODO: turn most of this into a dict
		self.aligned = -1

		self.qfullcov = -1
		self.sfullcov = -1
		self.fullcov = -1

		self.qtmcov = -1
		self.stmcov = -1
		self.tmcov = -1

		self.qlen = -1
		self.slen = -1

		self.qtmlen = -1
		self.stmlen = -1

		self.identity = -1

		self.distances = []
		self.qpresent = []
		self.spresent = []
		self.matrix = []

		self.quality = -1
		self.rmsd = -1
		self.length = -1

		self.name = ''

	def dump_json(self):
		self.IS_ANYTHING_WRONG()
		obj = {'query':self.query, 'subject': self.subject, 
			'queryfn':self.queryfn, 'subjectfn':self.subjectfn, 
			'qhel':self.qhel, 'shel':self.shel, 
			'qaligned':self.qaligned, 'saligned':self.saligned, 
			'qpresent':self.qpresent, 'spresent':self.spresent, 
			'distances':self.distances, 
			'matrix':self.matrix, 
			'quality':self.quality, 'rmsd':self.rmsd, 'length':self.length, 
			'aligned':self.aligned, 
			'qfullcov':self.qfullcov, 'sfullcov':self.sfullcov, 'fullcov':self.fullcov, 
			'qtmcov':self.qtmcov, 'stmcov':self.stmcov, 'tmcov':self.tmcov, 
			'qtmlen':self.qtmlen, 'stmlen':self.stmlen, 
			'qlen':self.qlen, 'slen':self.slen, 
			'identity':self.identity,
			'name':self.name,
		}
		return json.dumps(obj)

	def IS_ANYTHING_WRONG(self):
		assert self.query
		assert self.subject
		assert self.queryfn
		assert self.subjectfn
		assert self.matrix
		assert self.quality != -1
		assert self.rmsd != -1
		assert self.length != -1
		assert self.name
		return 0

	@staticmethod
	def parse_json(s):
		obj = json.loads(s)
		selfobj = Alignment()
		self.query = obj['query']
		self.subject = obj['subject']
		self.qhel = obj['qhel']
		selfobj.qaligned = obj['qaligned']
		selfobj.saligned = obj['saligned']
		selfobj.distances = obj['distances']
		selfobj.qpresent = obj['qpresent']
		selfobj.spresent = obj['spresent']
		selfobj.matrix = obj['matrix']
		selfobj.quality = obj['quality']
		selfobj.rmsd = obj['rmsd']
		selfobj.length = obj['length']
		selfobj.qlen = obj['qlen']
		selfobj.slen = obj['slen']
		selfobj.qcov = obj['qcov']
		selfobj.scov = obj['scov']
		return obj

	@staticmethod
	def parse_spout(spout):
		selfobj = Alignment()

		mode = 'fn'
		qlast = None
		slast = None
		for l in spout.split('\n'):
			if not l.strip(): continue
			elif l.lstrip().startswith('#'): continue
			elif mode == 'fn':
				if 'CID' in l:
					span = [int(i) for i in l.split()[-1][3:-1].split('-')]
					if not selfobj.qpresent: selfobj.qpresent.append(span)
					else: selfobj.spresent.append(span)
				elif l.startswith(' Query  '): 
					selfobj.queryfn = l[12:].strip()
				elif l.startswith(' and Target '): 
					selfobj.subjectfn = l[12:].strip()
					mode = 'matrix'
			elif mode == 'matrix':
				if l.strip().startswith('in fractional coordinates'): mode = 'scores'
				try: 
					row = [float(x) for x in l.split()]
					selfobj.matrix.append(row)
				except ValueError: pass
			elif mode == 'scores':
				if l.startswith('   quality'): selfobj.quality = float(l[15:22])
				elif l.startswith('     r.m.s.d'): selfobj.rmsd = float(l[15:22])
				elif l.startswith('      Nalign'): 
					selfobj.length = int(l[15:22])
					mode = 'prealignment'
			elif mode == 'prealignment':
				if l.startswith (' |---'): mode = 'alignment'
			elif mode == 'alignment':
				if l.startswith(' `----'): 
					mode = 'end'
					continue

				try: qcur = int(l[10:15].strip())
				except ValueError: qcur = None

				if selfobj.qaligned:
					if qlast is None and qcur is not None: selfobj.qaligned.append([qcur, qcur])
					elif qlast is None and qcur is None: selfobj.qaligned[-1][1] += 1
					elif qlast is not None and qcur is not None: 
						if qcur == qlast + 1: selfobj.qaligned[-1][1] += 1
						else: selfobj.qaligned.append([qcur, qcur])
					elif qlast is not None and qcur is None: selfobj.qaligned.append([qcur, 1])
				else:
					if qcur is None: selfobj.qaligned.append([qcur, 1])
					elif qcur is not None: selfobj.qaligned.append([qcur, qcur])
				qlast = qcur

				try: selfobj.distances.append(float(l[18:23].strip()))
				except ValueError: selfobj.distances.append(None)

				try: scur = int(l[35:39].strip())
				except ValueError: scur = None

				if selfobj.saligned:
					if slast is None and scur is not None: selfobj.saligned.append([scur, scur])
					elif slast is None and scur is None: selfobj.saligned[-1][1] += 1
					elif slast is not None and scur is not None: 
						if scur == slast + 1: selfobj.saligned[-1][1] += 1
						else: selfobj.saligned.append([scur, scur])
					elif slast is not None and scur is None: selfobj.saligned.append([scur, 1])
				else:
					if scur is None: selfobj.saligned.append([scur, 1])
					elif scur is not None: selfobj.saligned.append([scur, scur])
				slast = scur

		return selfobj
class Superpose(object):
	def __init__(self, d2dir='deuterocol2', force=False):
		self.d2dir = d2dir
		self.famdirs = []
		self.force = force
		if os.path.isdir('{}/config'.format(self.d2dir)) and os.path.isdir('{}/superpositions'.format(self.d2dir)):
			self.famdirs = [self.d2dir]
		else:
			for fn in os.listdir(self.d2dir):
				if fn.startswith('.'): continue
				elif os.path.isdir('{}/{}'.format(self.d2dir, fn)) and '_vs_' in fn:
					self.famdirs.append('{}/{}'.format(self.d2dir, fn))

	def main(self, famdir):
		done = []
		if VERBOSITY: info('Checking for existing alignments in {}...'.format(famdir))
		if not self.force and os.path.isfile('{}/superpositions/sp_all.tsv'.format(famdir)):
			with open('{}/superpositions/sp_all.tsv'.format(famdir)) as f:
				for l in f: 
					if not l.strip(): continue
					elif l.startswith('#'): continue
					else: done.append(l.split('\t')[0])
		if done: info('Skipping {} alignments (already done)'.format(len(done)))

		todo = -len(done)
		with open('{}/config/agenda.json'.format(famdir)) as g:
			for l in g: todo += 1

		n = 0
		t = time.time()
		if VERBOSITY: info('Performing {} alignments with SSM Superpose...'.format(todo))
		if todo:
			with open('{}/superpositions/sp_all.tsv'.format(famdir), 'a') as f:
				with open('{}/config/agenda.json'.format(famdir)) as g:
					for l in g:
						obj = json.loads(l)
						if not self.force and obj['name'] in done: continue
						if not n % 100: 
							info('Finished {}/{} alignments ({:0.2f}s since last message)'.format(n, todo, time.time() - t))
							t = time.time()
						cmd = ['superpose', \
							'{}/pdbs/{}.pdb'.format(famdir, obj['query']), \
							'-s', '{}/{}-{}'.format(obj['qchain'], *obj['qspan']), \
							'{}/pdbs/{}.pdb'.format(famdir, obj['subject']), \
							'-s', '{}/{}-{}'.format(obj['schain'], *obj['sspan']), \
						]
						out = subprocess.check_output(cmd)
						sp = Alignment.parse_spout(out)
						sp.query = obj['query']
						sp.queryfn = cmd[1]
						sp.qchain = obj['qchain']
						sp.qpresent = [obj['qspan']]
						sp.qhel = obj['qhelices']
						sp.subject = obj['subject']
						sp.subjectfn = cmd[4]
						sp.schain = obj['schain']
						sp.spresent = [obj['sspan']]
						sp.shel = obj['shelices']
						f.write('{}\t{}\n'.format(obj['name'], sp.dump_json()))
						n += 1
		info('Finished alignments for {}'.format(famdir))

	def run(self):
		if VERBOSITY: info('Checking for expected directory contents...')
		for famdir in self.famdirs:
			if not os.path.isfile('{}/config/agenda.json'.format(famdir)): raise IOError('Deuterocol2 directory structure not complete for {}'.format(famdir))

		for famdir in self.famdirs:
			lockfn = '{}/.lockfile'.format(famdir)
			skip = False
			try: 
				if os.path.isfile(lockfn): 
					with open(lockfn) as f:
						warn('Found lockfile in {} dated {}, skipping'.format(famdir, f.read().strip()))
						skip = True
						continue
				else:
					with open(lockfn, 'w') as f: f.write(time.strftime('%Y-%m-%d %H:%M:%S'))
				self.main(famdir)
			finally: 
				if os.path.isfile(lockfn) and not skip: os.remove(lockfn)
		info('Finished all assigned alignments')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('d2dir', default='deuterocol2', help='Deuterocol2 directory (contains config/, pdbs/, and superpositions/)')

	args = parser.parse_args()


	x = Superpose(d2dir=args.d2dir)
	x.run()
