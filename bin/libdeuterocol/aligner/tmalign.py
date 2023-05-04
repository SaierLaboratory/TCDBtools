import time
import sys
import os
import psutil
import shutil
import warnings
import subprocess
import json
import tempfile
import tqdm

def _parse_hrange(s): 
	try: return [int(x) for x in s[1:].split('-')]
	except ValueError: return None

def rle(arr):
	last = None
	newarr = []
	for elem in arr:
		if not newarr: newarr.append([elem, 1])
		elif last == elem: newarr[-1][-1] += 1
		else: newarr.append([elem, 1])
		last = elem
	return newarr

def rrle(arr):
	newarr = []
	for elem, count in arr:
		newarr.extend([elem] * count)
	return newarr

def compress_ranges(arr):
	last = None
	newarr = []
	for elem in arr:
		if not newarr: newarr.append([elem, 1])
		elif elem == last and elem is None: newarr[-1][-1] += 1
		elif elem is not None and last is None: newarr.append([elem, 1])
		elif elem == (last + 1): newarr[-1][-1] += 1
		else: newarr.append([elem, 1])
		last = elem
	return newarr

def decompress_ranges(arr):
	newarr = []
	for elem, count in arr:
		if elem is None: newarr.extend([None] * count)
		else: newarr.extend(range(elem, elem+count))
	return newarr

class TMAlign(object):
	def __init__(self, deut2, outdir):
		self.deut2 = deut2
		self.outdir = outdir

	def run(self):
		for pairbn in sorted(os.listdir(self.outdir)):
			pairdir = '{}/{}'.format(self.outdir, pairbn)
			if '_vs_' not in pairdir: continue
			if '.lockfile' in os.listdir(pairdir): 
				with open('{}/.lockfile'.format(pairdir)) as fh:
					pid, timestamp = fh.read().strip().split('\t')
					running = psutil.pid_exists(int(pid))
					if running: 
						warnings.warn('Found active lockfile. Skipping {}'.format(pairdir))
						continue
					else:
						warnings.warn('Found abandoned lockfile. Removing and continuing')
						os.remove('{}/.lockfile'.format(pairdir))

			try:
				with open('{}/.lockfile'.format(pairdir), 'w') as fh:
					out = '{pid}\t{datetime}\n'.format(pid=os.getpid(), datetime=time.strftime('%Y-%m-%d %H:%M:%S'))
					fh.write(out)

				if not os.path.isdir('{}/tmalignments'.format(pairdir)): os.mkdir('{}/tmalignments'.format(pairdir))

				total = len(list(self.deut2.get_alignments()))

				#TODO: Resume if halted early
				n = 0
				print('Performing alignments for {}'.format(pairdir))
				with open('{}/config/agenda.json'.format(pairdir)) as fh:
					with open('{}/tmalignments/sp_all.tsv'.format(pairdir), 'w') as outfh:
						localtotal = sum([1 for l in fh])
						fh.seek(0)
						for l in tqdm.tqdm(fh.read().split('\n')):
							if not l.strip(): continue
							#if not n % 100: print('Performed {}/{} alignments'.format(n, localtotal))
							inobj = json.loads(l)
							outobj = run_tmalign(inobj, pdbdir='{}/cut_pdbs'.format(self.outdir))
							qpdbc = '{query}_{qchain}'.format(**inobj)
							spdbc = '{subject}_{schain}'.format(**inobj)
							qhel = _parse_hrange(inobj['qlabel'][0])
							shel = _parse_hrange(inobj['slabel'][0])
							qlengths = self.deut2.get_pdblength(qpdbc, tmrange=qhel)
							slengths = self.deut2.get_pdblength(spdbc, tmrange=shel)
							outobj['qfullcov'] = outobj['length'] / qlengths['total']
							outobj['sfullcov'] = outobj['length'] / slengths['total']

							outobj['qtmaligned'] = self.count_aligned(self.deut2.get_states(qpdbc, tmrange=qhel), decompress_ranges(outobj['qpresent']), rrle(outobj['distances']))
							outobj['stmaligned'] = self.count_aligned(self.deut2.get_states(spdbc, tmrange=shel), decompress_ranges(outobj['spresent']), rrle(outobj['distances']))

							#FIXME: Loop covs
							outobj['qtmcov'] = outobj['qtmaligned'] / qlengths['tm']
							outobj['stmcov'] = outobj['stmaligned'] / slengths['tm']

							outfh.write('{}\t{}\n'.format(inobj['name'], json.dumps(outobj)))
							n += 1
			finally: os.remove('{}/.lockfile'.format(pairdir))

	def count_aligned(self, states, gappy, distances, threshold=5):
		return sum([state * (pos is not None) * (dist is not None and dist <= threshold) for state, pos, dist in zip(states, gappy, distances)])

	def tabulate(self):
		if not os.path.isdir('{}/tables'.format(self.outdir)): os.mkdir('{}/tables'.format(self.outdir))
		for pairbn in sorted(os.listdir(self.outdir)):
			pairdir = '{}/{}'.format(self.outdir, pairbn)
			if '_vs_' not in pairdir: continue
			with open('{}/tables/{}.tsv'.format(self.outdir, pairbn), 'w') as fh:
				self.tabulate_d2dir(pairdir, outfile=fh)

	def tabulate_d2dir(self, d2dir, outfile=None):
		if outfile is None: outfile = sys.stdout
		with open('{}/tmalignments/sp_all.tsv'.format(d2dir)) as fh:
			outfile.write('#name\tquery\tqchain\tqtms\tsubject\tschain\tstms\trmsd\tquality\tlength\tmaxtmcov\tmaxfullcov\n')
			for l in fh:
				name, jsonstr = l.split('\t')
				obj = json.loads(jsonstr)
				obj['name'] = name

				query, qchain, qtms, vs, subject, schain, stms = name.split('_')
				obj['query'] = query
				obj['qchain'] = qchain
				obj['qtms'] = qtms
				obj['subject'] = subject
				obj['schain'] = schain
				obj['stms'] = stms
				obj['maxtmcov'] = max(obj['qtmcov'], obj['stmcov'])
				obj['maxfullcov'] = max(obj['qfullcov'], obj['sfullcov'])

				outfile.write('{name}\t{query}\t{qchain}\t{qtms}\t{subject}\t{schain}\t{stms}\t{rmsd:0.2f}\t{quality:0.5f}\t{length}\t{maxtmcov:0.2%}\t{maxfullcov:0.2%}\n'.format(**obj))
				
def parse_matrix(fh):
	matrix = []
	for i, l in enumerate(fh): 
		if 2 <= i <= 4: matrix.append([float(x) for x in l.strip().split()[1:]])
	return matrix

def tmalign_crashed(e):
	obj = {}
	obj['qlength'] = 0
	obj['slength'] = 0
	obj['length'] = 0
	obj['rmsd'] = 999.0
	obj['quality'] = 0.0

	obj['distances'] = []
	obj['qpresent'] = []
	obj['spresent'] = []

	obj['matrix'] = [
		[1.0,0.0,0.0,0.0],
		[0.0,1.0,0.0,0.0],
		[0.0,0.0,1.0,0.0]
	]

	obj['error'] = str(e)
	return obj

def parse_tmalign(tmaout):
	obj = {}
	alnrows = []
	record = False
	obj['error'] = None

	for l in tmaout.split('\n'):
		if l.startswith('Length of Chain_1'): obj['qlength'] = int(l.split()[3])
		elif l.startswith('Length of Chain_2'): obj['slength'] = int(l.split()[3])
		elif l.startswith('Aligned length'):
			sl = l.replace(',' ,'').split()
			obj['length'] = int(sl[2])
			obj['rmsd'] = float(sl[4])

		elif l.startswith('TM-score='):
			tmscore = float(l.split()[1])
			if 'quality' in obj and tmscore > obj['quality']: obj['quality'] = tmscore
			else: obj['quality'] = tmscore

		elif l.startswith('(":"'): record = True
		elif record: alnrows.append(l.replace('\n', ''))

	for i, row in enumerate(alnrows[:3]):
		#FIXME: Everything is horrible
		if i == 0:
			qpresent = [j if row[j] != '-' else None for j in range(len(row))]
			obj['qpresent'] = compress_ranges(qpresent)
		elif i == 2:
			spresent = [j if row[j] != '-' else None for j in range(len(row))]
			obj['spresent'] = compress_ranges(spresent)
		elif i == 1:
			distances = []
			for c in row:
				if c == ' ': distances.append(None)
				elif c == ':': distances.append(0.0)
				else: distances.append(5.0)
			obj['distances'] = rle(distances)
			obj['truelength'] = len(row) - 1 - row.count(' ')

	return obj



def run_tmalign(obj, pdbdir):
	outobj = {}
	queryfn = '{pdbdir}/{query}_{qchain}_{qlabel}.pdb'.format(pdbdir=pdbdir, **obj)
	outobj['queryfn'] = queryfn
	subjectfn = '{pdbdir}/{subject}_{schain}_{slabel}.pdb'.format(pdbdir=pdbdir, **obj)
	outobj['subjectfn'] = subjectfn
	matrixfile = tempfile.NamedTemporaryFile()
	cmd = ['tmalign', queryfn, subjectfn, '-m', matrixfile.name]
	try: 
		out = subprocess.check_output(cmd)
		matrixfile.flush()
		matrixfile.seek(0)

		outobj.update(parse_tmalign(out.decode('utf-8')))
		outobj['matrix'] = parse_matrix(matrixfile)
	except subprocess.CalledProcessError as e:
		outobj.update(tmalign_crashed(e))

	return outobj

