from __future__ import division
import re
from Bio.Seq import Seq

class Alignment(object):
	opt = None
	zscore = None
	bitscore = None
	evalue = None
	score = None
	ident = None
	simil = None
	overlap = None
	qlen = None
	qoffset = None
	qstart = None
	qend = None
	tlen = None
	toffset = None
	tstart = None
	tend = None

	qseq = ''
	tseq = ''
	midline = ''

	def __init__(self):
		self.query = Seq('')
		#print(dir(self.query))
		self.target = Seq('')

	def emboss_style(self):
		qname = self.query.name[:13].ljust(13)
		tname = self.target.name[:13].ljust(13)
		out = ''

		out += '# \n'
		out += '# Aligned_sequences: 2\n'
		out += '# 1: {}\n'.format(self.query.name)
		out += '# 2: {}\n'.format(self.target.name)
		out += '# Matrix: {}\n'.format(self.matrix)
		out += '# Gap_penalty: {:0.1f}\n'.format(self.gapopen)
		out += '# Extend_penalty: {:0.1f}\n'.format(self.gapextend)
		out += '# \n'
		out += '# Length: {}\n'.format(self.overlap)
		out += '# Identity:'.ljust(17)
		out += '{}/{} ({:0.1%})\n'.format(int(self.ident * self.overlap), self.overlap, self.ident)
		out += '# Similarity:'.ljust(17)
		out += '{}/{} ({:0.1%})\n'.format(int(self.simil * self.overlap), self.overlap, self.simil)
		out += '# Gaps:'.ljust(17)
		out += '{}/{} ({:0.1%})\n'.format(self.midline.count('-'), self.overlap, self.midline.count('-')/self.overlap)
		out += '# Score: {:0.1f}\n'.format(self.score)
		out += '# \n'
		out += '# \n'
		out += '#{}\n'.format('='*39)

		out += '\n'

		def count_residues(s): return sum([r in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' for r in s])

		width = 50
		qi, ti = self.qstart, self.tstart
		for i in range(0, max(len(self.qseq), len(self.tseq), len(self.midline)), width):
			qseg = self.qseq[i:i+width]
			tseg = self.tseq[i:i+width]
			mseg = self.midline[i:i+width].replace(':', '|').replace('.', ':').replace('-', ' ')
			qdelta = count_residues(qseg)
			tdelta = count_residues(tseg)

			if len(qseg) < len(tseg): qseg += ' ' * (len(tseg) - len(qseg))
			elif len(tseg) > len(qseg): tseg += ' ' * (len(qseg) - len(tseg))

			out += qname + str(qi).rjust(7) + ' ' + qseg + str(qi+qdelta-1).rjust(7) + '\n'
			out += ' '*21 + mseg + '\n'
			out += tname + str(ti).rjust(7) + ' ' + tseg + str(ti+tdelta-1).rjust(7) + '\n'
			out += '\n'
	
			qi += qdelta
			ti += tdelta

		out += '\n'

		out += '#' + '-'*39


		return out

class AlignmentCollection(object):
	alignments = []
	def __init__(self):
		pass

	def __len__(self): return len(self.alignments)
	def __getitem__(self, index): return self.alignments.__getitem__(index)
	def __iter__(self): return iter(self.alignments)

	def get_opt(self): return [aln.opt for aln in self]
	def get_zscore(self): return [aln.zscore for aln in self]
	def get_bitscore(self): return [aln.bitscore for aln in self]
	def get_evalue(self): return [aln.evalue for aln in self]
	def get_score(self): return [aln.score for aln in self]
	def get_ident(self): return [aln.ident for aln in self]
	def get_simil(self): return [aln.simil for aln in self]
	def get_overlap(self): return [aln.overlap for aln in self]
	def get_length(self): return [aln.overlap for aln in self]
	def get_midline(self): return [aln.midline for aln in self]

	@staticmethod
	def parse(f, mode=10):

		obj = AlignmentCollection()

		matrix = None
		gapopen = None
		gapextend = None

		if mode == 10:

			def _split_dataline(l):
				if l.startswith(';'): return re.split('\s*:\s+', l[2:].strip())
				else: raise ValueError('Invalid dataline (Is this really a "-m 10" file?)')

			state = ''
			alignment = None
			for l in f:

				if l.startswith('>>>'): state = 'prog'

				elif l.startswith('>>'): 
					if alignment is not None: obj.alignments.append(alignment)
					state = 'alignment'
					alignment = Alignment()
					alignment.target.name = l[2:].strip()
					alignment.matrix = matrix
					alignment.gapopen = gapopen
					alignment.gapextend = gapextend

				elif l.startswith('>'): 
					if state == 'alignment': 
						state = 'query'
						alignment.query.name = l[1:].strip()
					elif state == 'query': 
						state = 'target'
						#alignment.target.name = l[1:]

				else:
					if state == 'prog':
						if 'pg_matrix' in l:
							k, v = _split_dataline(l)
							matrix = v
						elif 'pg_open-ext' in l:
							k, v = _split_dataline(l)
							gapopen, gapextend = [-int(x) for x in v.split()]

					elif state == 'alignment':
						k, v = _split_dataline(l)
						if k.endswith('opt'): alignment.opt = int(v)
						elif k.endswith('z-score'): alignment.zscore = float(v)
						elif k.endswith('bits'): alignment.bitscore = float(v)
						elif k.endswith('expect'): alignment.evalue = float(v)
						elif k.endswith('score'): alignment.score = int(v)
						elif k.endswith('ident'): alignment.ident = float(v)
						elif k.endswith('sim'): alignment.simil = float(v)
						elif k.endswith('overlap'): alignment.overlap = int(v)

					elif state == 'query':
						try: 
							k, v = _split_dataline(l)
							if k.endswith('len'): alignment.qlen = int(v)
							elif k.endswith('offset'): alignment.qoffset = int(v)
							elif k.endswith('al_start'): alignment.qstart = int(v)
							elif k.endswith('al_stop'): alignment.qend = int(v)
						except ValueError:
							alignment.qseq += l.replace('\n', '')

					elif state == 'target':
						if 'al_cons' in l: state = 'midline'
						else:
							try: 
								k, v = _split_dataline(l)
								if k.endswith('len'): alignment.tlen = int(v)
								elif k.endswith('offset'): alignment.toffset = int(v)
								elif k.endswith('start'): alignment.tstart = int(v)
								elif k.endswith('stop'): alignment.tend = int(v)
							except ValueError:
								alignment.tseq += l.replace('\n', '')

					elif state == 'midline':
						alignment.midline += l.replace('\n', '')

			if alignment is not None: obj.alignments.append(alignment)

			return AlignmentCollection()
		else: raise NotImplementedError('Mode {} parser not implemented (use -m 10)'.format(mode))
			
