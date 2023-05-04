#!/usr/bin/env python2

from __future__ import print_function, division, generators

import argparse, json, os, subprocess, re, sys, time
import tempfile
import superpose
import shutil

import tmalignparser


VERBOSITY = 1
def info(*things):
	print('[INFO]:', *things, file=sys.stderr)
def warn(*things):
	print('[WARNING]:', *things, file=sys.stderr)
def error(*things):
	print('[ERROR]:', *things, file=sys.stderr)
	exit(1)

def rle(arr):
	newarr = []
	lastx = 'INVALID'
	for x in arr:
		if lastx == 'INVALID': newarr.append([x, 1])
		else:
			if x == newarr[-1][0]: newarr[-1][1] += 1
			else: newarr.append([x, 1])
		lastx = x
	return newarr

def count_residues(pdbfn):
	n = 0
	with open(pdbfn) as f:
		for l in f:
			if l.startswith('ATOM'): 
				if l[13:15] == 'CA': n += 1
	return n

def manually_cut(ifn, ofn, chain, start, end):

	intended = (start, end)
	def intersects(a, b):
		if (a[0] < b[0] and a[1] < b[0]): return False
		elif (a[0] > b[1] and a[1] > b[1]): return False
		else: return True

	#skip negative residues for now because TMalign crashes on them
	start = max(0, start)
	
	with open(ifn) as infile:
		with open(ofn, 'w') as outfile:
			for l in infile:
				if l.startswith('DBREF '):
					if l[12] != chain: continue
					elif l[14:18].strip(): 
						thisref = (int(l[14:18]), int(l[20:24]))
						if not intersects(thisref, intended): continue
				elif l.startswith('SEQRES'):
					if l[11] != chain: continue
				elif l.startswith('HET   '):
					if l[12] != chain: continue
				#elif l.startswith('ANISOU'): continue
				elif l.startswith('ATOM  '):
					if l[21] != chain: continue
					resi = int(l[22:26])
					if not (start <= resi <= end): continue
				elif l.startswith('TER   '):
					if l[21] != chain: continue
					resi = int(l[22:26])
					if not (start <= resi <= end): continue
				elif l.startswith('HETATM'):
					if l[21] != chain: continue
					resi = int(l[22:26])
					if not (start <= resi <= end): continue
				elif l.startswith('ENDMDL'):
					outfile.write(l)
					outfile.write('END' + 77*' ' + '\n')
					break
				outfile.write(l)

def truncate_to_resolved(pdbfn, chain, start, end):
	resolved = []
	with open(pdbfn) as f:
		for l in f:
			if l.startswith('ATOM'):
				if l[13:15] == 'CA':
					if (l[21] == chain) or ((l[21] == ' ') and (chain == 'A')):
						resolved.append(int(l[22:26]))
	newstart = None
	newend = None
	for i in range(start, end+1):
		if i in resolved: 
			if newstart is None: newstart = i
			else: newstart = min(newstart, i)
			if newend is None: newend = i
			else: newend = max(newend, i)

	return newstart - start, end - newend

def _intersection_size(spans1, spans2):
	n = 0
	lastn = None

	#TODO: optimize to switch spans1 to shorterspans
	for span1 in spans1:
		ind1 = set(range(span1[0], span1[1]+1))
		for span2 in spans2:
			ind2 = set(range(span2[0], span2[1]+1))
			n += len(ind1.intersection(ind2))
	return n


class TMalign(superpose.Superpose):
	def __init__(self, d2dir='deuterocol2', force=False, skip_cut=True, min_tms=3, ignore_empty=False):
		superpose.Superpose.__init__(self, d2dir=d2dir, force=force)
		self.skip_cut = skip_cut
		self.min_tms = min_tms
		self.ignore_empty = ignore_empty

	def main(self, famdir):
		done = []
		if not os.path.isdir('{}/tmalignments'.format(famdir)): os.mkdir('{}/tmalignments'.format(famdir))
		if VERBOSITY: info('Checking for existing alignments in {}...'.format(famdir))
		if not self.force and os.path.isfile('{}/tmalignments/sp_all.tsv'.format(famdir)):
			with open('{}/tmalignments/sp_all.tsv'.format(famdir)) as f:
				for l in f:
					if not l.strip(): continue
					elif l.startswith('#'): continue
					else: 
						try: 
							json.loads(l.split('\t')[1])
							done.append(l.split('\t')[0])
						except ValueError: break
		if done: info('Skipping {} alignments (already done)'.format(len(done)))
			

		todo = -len(done)
		with open('{}/config/agenda.json'.format(famdir)) as g:
			for l in g: 
				if not l.strip(): continue
				elif l.startswith('#'): continue
				todo += 1
		n = 0
		t = superpose.time.time()
		logfile = open('{}/errors.log'.format(famdir), 'w')
		if not todo:
			with open('{}/DONE'.format(famdir), 'w') as fh: pass
		else:
			if VERBOSITY: info('Performing {} alignments with TM-align...'.format(todo))
			with open('{}/tmalignments/sp_all.tsv'.format(famdir), 'a') as f:
				with open('{}/config/agenda.json'.format(famdir)) as g:
					for l in g:
						alnstart = time.time()
						obj = json.loads(l)
						if not self.force and obj['name'] in done: 
							n += 1
							continue
						if not n % 100:
							#TODO: fix this progress bar thing
							info('Finished {}/{} alignments ({:0.2f}s since last message)'.format(n, todo, superpose.time.time() - t))
							t = superpose.time.time()
						if len(obj['qindices']) < self.min_tms: 
							n += 1
							continue
						if len(obj['sindices']) < self.min_tms: 
							n += 1
							continue
						query, subject = obj['name'].split('_vs_')
						tf = tempfile.NamedTemporaryFile()
						qfn = '{}/../cut_pdbs/{}.pdb'.format(famdir, query)
						sfn = '{}/../cut_pdbs/{}.pdb'.format(famdir, subject)
						cmd = ['TMalign', 
							qfn, sfn, 
							'-m', tf.name
						]

						#this is where TMalign is run
						p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
						out, err = p.communicate()
						if p.returncode:
							#29: file not found, e.g. from being incompetently eliminated from the pipeline
							#174: segfault, seems to be from having empty structures
							#if p.returncode in {29, 174}: 
							if p.returncode in {174}:
								n += 1
								if self.ignore_empty: continue
								warn('{} or {} caused a problem!'.format(qfn, sfn))
								logfile.write('{} or {} caused a problem (174)!\n'.format(qfn, sfn))
								continue
							elif p.returncode in {29}:
								warn('{} or {} does not exist'.format(qfn, sfn))
								logfile.write('{} or {} does not exist (29)\n'.format(qfn, sfn))
								continue
							else: 
								print('CMD')
								print(' '.join(cmd))
								print('STDOUT')
								print(out)
								print('STDERR')
								print(err)
								print('Error: {name}: {err}, return code = {returncode}'.format(name=obj['name'], err=err, returncode=p.returncode))
								logfile.write('{} vs. {} had an unknown problem (-11)\n'.format(qfn, sfn))
								continue
								#raise Exception('Error: {err}, return code = {returncode}'.format(err=err, returncode=p.returncode))
						tf.flush()
						tf.seek(0)

						#capture stderr to do this later #TODO
						#if out.startswith('forrtl: No such file'): 
						#	superpose.error('Could not find at least one of {} and {}'.format(qfn, tfn))
						#for l in out.split('\n'): print(l)
						sp = superpose.Alignment()
						for i, rawrow in enumerate(tf):
							if re.match('^ ?[0-3]', rawrow):
								if len(sp.matrix) < 3:
									row = [float(x) for x in rawrow.strip().split()[1:]]
									sp.matrix.append(row[1:] + row[:1])
						sp.qpresent = [obj['qspan']]
						sp.spresent = [obj['sspan']]
						sp.query = obj['query']
						sp.subject = obj['subject']
						sp.queryfn = qfn
						sp.subjectfn = sfn
						#print(sp.dump_json())
						tmscores = []
						alignment = False
						aligned = []

						#this is the where TMalign parsing starts
						for tml in out.split('\n'):
							if alignment is False:
								if tml.startswith('Aligned length'):
									sp.length = int(tml.split()[2][:-1])
									sp.rmsd = float(tml.split()[4][:-1])
									sp.identity = float(tml.strip().split()[-1])
								elif tml.startswith('TM-score'):
									tmscores.append(float(tml.split()[1]))
								elif tml.startswith('(":" denotes'): alignment = True
							else:
								if 'Total running time' in tml: continue
								elif tml.strip():
									aligned.append(tml.replace('\n', ''))

						qseq, dseq, sseq = aligned
						sp.quality = max(tmscores)

						covstuff = True
						if covstuff:
							covstart = time.time()

							#if a contig couldn't be found but the fragment was shorter than this, skip
							minl_relevant = 4

							qatomseq = tmalignparser.extract_atom_sequence(open(qfn), obj['qchain'])
							satomseq = tmalignparser.extract_atom_sequence(open(sfn), obj['schain'])

							permitted = '.:'
							lastd = '#'
							qfrags = []
							sfrags = []
							for q, d, s in zip(qseq, dseq, sseq):
								if d in permitted: 
									if lastd in permitted:
										qfrags[-1] += q
										sfrags[-1] += s
									else:
										qfrags.append(q)
										sfrags.append(s)
								lastd = d

							qistart = None
							qaln_combined = tmalignparser.NumberedSequence()
							qcontigs = []
							for qfrag in qfrags:
								qcontig = qatomseq.gapless_align_string(qfrag, start=qistart)
								qcontigs.append(qcontig)
								if qcontig is None: 
									#if len(qfrag) < minl_relevant: continue
									#print(qfrag, qatomseq)
									continue

								qaln_combined = qaln_combined + qcontig
								qistart = qcontig.get_range()[-1] + 1
							#print(qcontigs)
							#print(qaln_combined)
							#print([len(x) for x in qcontigs], '=', len(qaln_combined))

							sistart = None
							saln_combined = tmalignparser.NumberedSequence()
							for sfrag in sfrags:
								scontig = satomseq.gapless_align_string(sfrag, start=sistart)
								#print(sfrag, sistart, scontig, satomseq)
								if scontig is None:
									#if len(sfrag) < minl_relevant: continue
									continue

								saln_combined = saln_combined + scontig
								#sistart = scontig.get_range()[-1] + 1
								sistart = scontig.get_range()[-1] + 1

							try:
								qindices = tmalignparser.collapse(obj['qindices'])
								sindices = tmalignparser.collapse(obj['sindices'])
							except IndexError:
								print(obj['qindices'])
								print(obj['sindices'])
								exit(1)
							#print(qaln_combined.get_ranges())
							#print(tmalignparser.collapse(qaln_combined.get_ranges()))
							n_qaligned = _intersection_size(tmalignparser.collapse(qaln_combined.get_ranges()), qindices)
							n_saligned = _intersection_size(tmalignparser.collapse(saln_combined.get_ranges()), sindices)
							#print(tmalignparser.collapse(saln_combined.get_ranges()), tmalignparser.collapse(obj['sindices']))
							#print(obj['slen'])
							#print(obj)

							sp.qlen = len(qatomseq)
							sp.slen = len(satomseq)

							#sp.qtmlen = obj['qlen']
							sp.qtmlen = _intersection_size(obj['qindices'], qindices)
							#print('q', n_qaligned, sp.qtmlen)
							#print('q', n_qaligned, sp.qtmlen)
							#sp.stmlen = obj['slen']
							sp.stmlen = _intersection_size(obj['sindices'], sindices)
							#print('s', n_saligned, sp.stmlen)
							#print('s', n_saligned, sp.stmlen)

							if n_qaligned > sp.qtmlen: 
								print(n_qaligned, sp.qtmlen)
								print(obj)
								raise IndexError('whyyyyyyyyyyyyyyyyyyyyyyyyyy')
							if n_saligned > sp.stmlen: 
								print(n_saligned, sp.stmlen)
								print(obj)
								raise IndexError('whyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy')

							sp.qtmcov = n_qaligned / sp.qtmlen
							sp.stmcov = n_saligned / sp.stmlen
							sp.tmcov = max(sp.qtmcov, sp.stmcov)

							sp.qfullcov = sp.length / sp.qlen
							sp.sfullcov = sp.length / sp.slen
							sp.fullcov = max(sp.qfullcov, sp.sfullcov)


							covtime = time.time() - covstart
						dseq += ' ' * max(0, len(qseq)-len(dseq))
						qaligned, saligned, distances = [], [], []
						aligned = 0


						lastqr, lastsr = '', ''
						lastdist = -1
						qi, si = sp.qpresent[0][0], sp.spresent[0][0]
						for qcur, dist, scur in zip(qseq, dseq, sseq):
							if lastdist  == -1:
								if dist == ' ': 
									qaligned.append([None, 1])
									saligned.append([None, 1])
									distances.append(None)
								elif dist == '.':
									if qcur == '-': qaligned.append([None, 1])
									else: qaligned.append([qi, 1])
									if scur == '-': saligned.append([None, 1])
									else: saligned.append([si, 1])
									distances.append(5.0)
									aligned += 1
								elif dist == ':':
									if qcur == '-': qaligned.append([None, 1])
									else: qaligned.append([qi, 1])
									if scur == '-': saligned.append([None, 1])
									else: saligned.append([si, 1])
									distances.append(0.0)
									aligned += 1
							elif dist == ' ': 
								if qaligned[-1][0] is None: qaligned[-1][1] += 1
								else: qaligned.append([None, 1])
								if saligned[-1][0] is None: saligned[-1][1] += 1
								else: saligned.append([None, 1])
								distances.append(None)
							elif dist == '.': 
								if qaligned[-1][0] is None: qaligned.append([qi, 1])
								else: qaligned[-1][1] += 1
								if saligned[-1][0] is None: saligned.append([si, 1])
								else: saligned[-1][1] += 1
								distances.append(5.0)
								aligned += 1
							elif dist == ':': 
								if qaligned[-1][0] is None: qaligned.append([qi, 1])
								else: qaligned[-1][1] += 1
								if saligned[-1][0] is None: saligned.append([si, 1])
								else: saligned[-1][1] += 1
								distances.append(0.)
								aligned += 1
							#if not qaligned and qcur == '-': qaligned.append([None, 1])
							#elif not qaligned and qcur != '-': qaligned.append([qi, 1])
							#elif qcur == '-' and lastqr == '-': qaligned[-1][1] += 1
							#elif qcur == '-' and lastqr != '-': qaligned.append([None, 1])
							#elif qcur != '-' and lastqr != '-': qaligned[-1][1] += 1
							#elif qcur != '-' and lastqr == '-': qaligned.append([qi, 1])

							#if not saligned and scur == '-': saligned.append([None, 1])
							#elif not saligned and scur != '-': saligned.append([si, 1])
							#elif scur == '-' and lastsr == '-': saligned[-1][1] += 1
							#elif scur == '-' and lastsr != '-': saligned.append([None, 1])
							#elif scur != '-' and lastsr != '-': saligned[-1][1] += 1
							#elif scur != '-' and lastsr == '-': saligned.append([si, 1])
						
							##TODO: expose distance cutoff
							#if dist == ' ': distances.append(None)
							#elif dist == '.': distances.append(5.0)
							#elif dist == ':': distances.append(0.)

							if qcur != '-': qi += 1
							if scur != '-': si += 1
							if qcur == '-': qi += 1
							if scur == '-': si += 1
							lastqr = qcur
							lastsr = scur
							lastdist = dist

						sp.qaligned = qaligned
						sp.saligned = saligned
						sp.distances = rle(distances)
						sp.aligned = aligned

						query, qchain, qhel, vs, subject, schain, shel = obj['name'].split('_')
						sp.qhel = [int(x) for x in qhel[1:].split('-')]
						sp.shel = [int(x) for x in shel[1:].split('-')]

						sp.name = obj['name']

						f.write('{}\t{}\n'.format(obj['name'], sp.dump_json()))
						n += 1
						alntime = time.time() - alnstart
						#print(covtime, alntime)
		info('Finished alignments for {}'.format(famdir))


	def cut_pdbs(self, famdir):
		indices = {}
		with open('{}/config/indices.json'.format(famdir)) as f:
			for l in f:
				sl = l.split('\t')
				indices[sl[0]] = json.loads(sl[1])

		bundle = -1
		with open('{}/config/agenda.json'.format(famdir)) as f:
			obj = json.loads(f.readline())
			bundle = obj['bundle']

		cutme = set()
		with open('{}/config/align_me.json'.format(famdir)) as f:
			for l in f:
				obj = json.loads(l)
				for fam in obj:
					for pdbid in obj[fam]:
						cutme.add(pdbid[:4])

		chains = {}
		for pdbc in indices:
			try: chains[pdbc[:4]].append(pdbc)
			except KeyError: chains[pdbc[:4]] = [pdbc]

		with open('{}/../deuterocol1/termini.json'.format(famdir)) as fh:
			termini = json.load(fh)
		
		for pdb in sorted(cutme):
			for pdbc in chains[pdb]:
				#if len(indices[pdbc]) <= 2: continue
				#if bundle >= len(indices[pdbc]): 
				if False:
					shutil.copy('{}/../pdbs/{}.pdb'.format(famdir, pdb), 
						'{}/../cut_pdbs/{}_h{}-{}.pdb'.format(famdir, pdbc, 1, len(indices[pdbc])))
				elif len(indices[pdbc]) == 0 or bundle == 0:
					start, end = termini[pdbc[:4]][pdbc[-1]]
					infile = '{}/../pdbs/{}.pdb'.format(famdir, pdb)
					outfile = '{}/../cut_pdbs/{}_h1-1.pdb'.format(famdir, pdbc)
					if not os.path.isfile(infile): raise IOError('Could not find {}'.format(infile))
					if os.path.isfile(outfile) and os.path.getsize(outfile): continue
					manually_cut(infile, outfile, pdbc[-1], start, end)
				else:
					if len(indices[pdbc]) <= bundle:
						infile = '{}/../pdbs/{}.pdb'.format(famdir, pdb)
						outfile = '{}/../cut_pdbs/{}_h{}-{}.pdb'.format(famdir, pdbc, 1, len(indices[pdbc]))
						if not os.path.isfile(infile): raise IOError('Could not find {}'.format(infile))
						if os.path.isfile(outfile) and os.path.getsize(outfile): continue
						#print('pdbc={}, indices={}'.format(pdbc, indices[pdbc]))
						manually_cut(infile, outfile, chain=pdbc[-1], start=indices[pdbc][0][0], end=indices[pdbc][-1][1])
					else:
						for start in range(0, len(indices[pdbc])-bundle+1):
							end = start + bundle - 1
							infile = '{}/../pdbs/{}.pdb'.format(famdir, pdb)
							outfile = '{}/../cut_pdbs/{}_h{}-{}.pdb'.format(famdir, pdbc, start+1, end+1)
							if not os.path.isfile(infile): raise IOError('Could not find {}'.format(infile))
							if os.path.isfile(outfile) and os.path.getsize(outfile): continue

							manually_cut(infile, outfile, chain=pdbc[-1], start=indices[pdbc][start][0], end=indices[pdbc][end][1])


	def run(self):
		if VERBOSITY: info('Checking for expected directory contents...')
		for famdir in self.famdirs:
			if not os.path.isfile('{}/config/agenda.json'.format(famdir)): raise IOError('Missing agenda.json: Deuterocol2 directory structure not complete for {}'.format(famdir))
			if not os.path.isfile('{}/config/indices.json'.format(famdir)): raise IOError('Missing indices.json: Deuterocol2 directory structure not complete for {}'.format(famdir))

			if not os.path.isdir('{}/../cut_pdbs'.format(famdir)):
				os.mkdir('{}/../cut_pdbs'.format(famdir))
			if not self.skip_cut: 
				if VERBOSITY: info('Cutting PDBs for {}'.format(famdir))
				self.cut_pdbs(famdir)

		for famdir in self.famdirs:
			lockfn = '{}/.lockfile'.format(famdir)
			skip = False
			try: 
				if os.path.isfile('{}/DONE'.format(famdir)): 
					info('{} is done, skipping'.format(famdir))
					continue

				if os.path.isfile(lockfn): 
					with open(lockfn) as f:
						warn('Found lockfile in {} dated {}, skipping'.format(famdir, f.read().strip()))
						skip = True
						continue
				else:
					try:
						with open(lockfn, 'w') as f: f.write(time.strftime('%Y-%m-%d %H:%M:%S'))
					except IOError:
						warn('Could not find {}, skipping to next library pair'.format(famdir))
				self.main(famdir)
				with open('{}/DONE'.format(famdir), 'w') as fh: pass
			finally: 
				if os.path.isfile(lockfn) and not skip: os.remove(lockfn)
		info('Finished all assigned alignments')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--skip-cut', action='store_true')
	parser.add_argument('--min-tms', default=3, type=int, help='Minimum TMSs allowed in an aligned PDB')
	parser.add_argument('--ignore-empty', action='store_true', help='Ignore TMalign segfaults (mostly caused by ATOMless structures)')
	parser.add_argument('d2dir', default='deuterocol2', help='Deuterocol2 directory (contains config/, pdbs/, and superpositions/)')

	args = parser.parse_args()


	x = TMalign(d2dir=args.d2dir, skip_cut=args.skip_cut, min_tms=args.min_tms, ignore_empty=args.ignore_empty)
	x.run()
