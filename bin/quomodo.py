#!/usr/bin/env python

import gzip
import base64
import argparse
import tempfile
import threading
import time
import sys
import os
import io
import json
import subprocess
import re

try: import pymol
except ModuleNotFoundError: pymol = None
except ImportError: pymol = None

def warn(*things):
	print('[WARNING]', *things, file=sys.stderr)

def error(*things):
	print('[ERROR]', *things, file=sys.stderr)
	exit(1)

class Alignment(object):
	data = None
	deutdir = None

	def get_query_fn(self, cut=False): 
		if cut: return '{}/cut_pdbs/{}'.format(self.deutdir, os.path.basename(self.data['queryfn']))
		else: return '{}/pdbs/{}.pdb'.format(self.deutdir, os.path.basename(self.data['query']))
	def get_sbjct_fn(self, cut=False): 
		if cut: return '{}/cut_pdbs/{}'.format(self.deutdir, os.path.basename(self.data['subjectfn']))
		else: return '{}/pdbs/{}.pdb'.format(self.deutdir, os.path.basename(self.data['subject']))

	@staticmethod
	def find(deutdir, name):
		bn = os.path.basename(deutdir)
		found = False
		if '_vs_' in bn:
			actualdeutdir = '{}/..'.format(deutdir)
			if not os.path.isfile('{}/tmalignments/sp_all.tsv'.format(deutdir)): 
				error('No results recorded for {}'.format(bn))
			with open('{}/tmalignments/sp_all.tsv'.format(deutdir)) as fh:
				for l in fh:
					if l.startswith(name + '\t'): 
						found = l
						fampair = bn
						break
		else:
			actualdeutdir = deutdir
			for subdir in os.listdir(deutdir):
				if '_vs_' in subdir:
					if not os.path.isfile('{}/{}/tmalignments/sp_all.tsv'.format(deutdir, subdir)): 
						warn('No results recorded for {}'.format(subdir))
						continue
					with open('{}/{}/tmalignments/sp_all.tsv'.format(deutdir, subdir)) as fh:
						for l in fh:
							if l.startswith(name + '\t'): 
								found = l
								fampair = subdir
								break
		if not found: error('Could not find {} in {}'.format(name, deutdir))

		sp = Alignment()
		sp.data = json.loads(found[found.find('\t')+1:])
		sp.deutdir = actualdeutdir
		sp.fampair = fampair

		return sp

	def get_aligned(self, tmalign='TMalign'):
		query = self.get_query_fn(cut=True)
		sbjct = self.get_sbjct_fn(cut=True)

		matrixfile = tempfile.NamedTemporaryFile(mode='r')
		out = subprocess.check_output([tmalign, query, sbjct, '-m', matrixfile.name]).decode('utf-8')
		matrixfile.flush()
		matrixfile.seek(0)


		matrix = []
		for l in matrixfile:
			if re.match('^ ?[0-3]', l):
				if len(matrix) < 3:
					row = [float(x) for x in l.strip().split()[1:]]
					matrix.append(row[1:] + row[:1])

		qmasked = []
		aligned = []
		smasked = []
		recordthis = 0
		for l in out.split('\n'): 
			if 'denotes aligned residue' in l:
				recordthis = 3
			elif recordthis > 0:
				if not qmasked:
					j = 0
					for i, r in enumerate(l):
						if r == '-': qmasked.append(None)
						else: 
							qmasked.append(j)
							j += 1
				elif not aligned:
					for i, r in enumerate(l):
						if r == ' ': aligned.append(None)
						else: aligned.append(i)
				elif not smasked:
					j = 0
					for i, r in enumerate(l):
						if r == '-': smasked.append(None)
						else: 
							smasked.append(j)
							j += 1
				recordthis -= 1

		with open(query) as fh: qresilist = get_resilist(fh.read())
		with open(sbjct) as fh: sresilist = get_resilist(fh.read())


		qaligned = []
		saligned = []
		for i, r in enumerate(aligned):
			if r is None: pass
			else:
				qaligned.append(qresilist[qmasked[r]])
				saligned.append(sresilist[smasked[r]])

		return rlencode(qaligned), rlencode(saligned), matrix

def zipstr(s):
	fo = io.BytesIO()
	gzo = gzip.GzipFile(fileobj=fo, mode='wb')
	gzo.write(s.encode('utf-8'))
	gzo.flush()
	gzo.close()
	fo.flush()
	fo.seek(0)
	return fo.read()

def unzipstr(b):
	fo = io.BytesIO(b)
	gzo = gzip.GzipFile(fileobj=fo, mode='rb')
	return gzo.read()

def main_local(args): 
	deutdir = args.local[0]
	name = args.local[1]
	if pymol is None: raise ImportError('No module named \'pymol\'; Check PYTHONPATH?')
	alignment = Alignment.find(deutdir, name)
	qaligned, saligned, matrix = alignment.get_aligned(tmalign=args.tmalign)
	with open(alignment.get_query_fn()) as fh: querypdb = fh.read()
	with open(alignment.get_sbjct_fn()) as fh: sbjctpdb = fh.read()
	with open(alignment.get_query_fn(cut=True)) as fh: cutquerypdb = fh.read()
	with open(alignment.get_sbjct_fn(cut=True)) as fh: cutsbjctpdb = fh.read()
	data['querypdb'] = querypdb
	data['sbjctpdb'] = sbjctpdb
	data['cutquerypdb'] = cutquerypdb
	data['cutsbjctpdb'] = cutsbjctpdb
	data['alignment'] = alignment.data
	data['alignment']['qaligned'] = qaligned
	data['alignment']['saligned'] = saligned
	data['alignment']['matrix'] = matrix
	data['fampair'] = alignment.fampair

	dsthread = DoStuffThread(kwargs=dict(
		obj=obj,
		delay=args.delay,
	))
	dsthread.start()
	pymol.finish_launching()

	data = {}

def main_server(args):
	deutdir = args.server[0]
	name = args.server[1]
	alignment = Alignment.find(deutdir, name)
	qaligned, saligned, matrix = alignment.get_aligned(tmalign=args.tmalign)

	data = {}
	with open(alignment.get_query_fn()) as fh: querypdb = fh.read()
	with open(alignment.get_sbjct_fn()) as fh: sbjctpdb = fh.read()
	with open(alignment.get_query_fn(cut=True)) as fh: cutquerypdb = fh.read()
	with open(alignment.get_sbjct_fn(cut=True)) as fh: cutsbjctpdb = fh.read()
	data['querypdb'] = querypdb
	data['sbjctpdb'] = sbjctpdb
	data['cutquerypdb'] = cutquerypdb
	data['cutsbjctpdb'] = cutsbjctpdb
	data['alignment'] = alignment.data
	data['alignment']['qaligned'] = qaligned
	data['alignment']['saligned'] = saligned
	data['alignment']['matrix'] = matrix
	data['fampair'] = alignment.fampair
	jsonstr = json.dumps(data)
	zipped = zipstr(jsonstr)
	encoded = base64.b64encode(zipped).decode('utf-8')

	print('+++BEGIN+++')
	print(encoded)
	print('+++END+++')
	
def main_client(args):
	if pymol is None: raise ImportError('No module named \'pymol\'; Check PYTHONPATH?')

	relevant = False
	payload = ''
	for l in sys.stdin:
		if not relevant:
			if l.startswith('+++BEGIN+++'): relevant = True
		else:
			if l.startswith('+++END+++'): relevant = False
			else: payload += l
	decoded = base64.b64decode(payload)
	unzipped = unzipstr(decoded)
	obj = json.loads(unzipped)

	dsthread = DoStuffThread(kwargs=dict(
		obj=obj,
		delay=args.delay,
	))
	dsthread.start()
	pymol.finish_launching()

class DoStuffThread(threading.Thread):
	def run(self):
		time.sleep(self._kwargs['delay'])
		obj = self._kwargs['obj']
		fampair = obj['fampair']
		queryfh = tempfile.NamedTemporaryFile(suffix='.pdb')
		queryfh.write(obj['querypdb'].encode('utf-8'))
		queryfh.flush()
		cutqueryfh = tempfile.NamedTemporaryFile(suffix='.pdb')
		cutqueryfh.write(obj['cutquerypdb'].encode('utf-8'))
		cutqueryfh.flush()
		queryname = os.path.basename(obj['alignment']['queryfn'])
		queryname = queryname[:queryname.find('.')]

		sbjctfh = tempfile.NamedTemporaryFile(suffix='.pdb')
		sbjctfh.write(obj['sbjctpdb'].encode('utf-8'))
		sbjctfh.flush()
		#cutsbjctfh = tempfile.NamedTemporaryFile(suffix='.pdb')
		#cutsbjctfh.write(obj['cutsbjctpdb'].encode('utf-8'))
		#cutsbjctfh.flush()
		sbjctname = os.path.basename(obj['alignment']['subjectfn'])
		sbjctname = sbjctname[:sbjctname.find('.')]

		tmpqueryname = os.path.basename(queryfh.name).replace('.pdb', '')
		tmpcutqueryname = os.path.basename(cutqueryfh.name).replace('.pdb', '')
		tmpsbjctname = os.path.basename(sbjctfh.name).replace('.pdb', '')
		#tmpcutsbjctname = os.path.basename(cutsbjctfh.name).replace('.pdb', '')

		matrix = obj['alignment']['matrix']
		pymolmatrix = []
		for row in matrix: pymolmatrix.extend(row)
		pymolmatrix.extend([0,0,0,1])

		pymol.cmd.load(cutqueryfh.name)
		pymol.cmd.load(queryfh.name)

		#pymol.cmd.load(cutsbjctfh.name)
		pymol.cmd.load(sbjctfh.name)

	
		print(summary_stats(obj, name='{}/{}'.format(fampair, obj['alignment'].get('name', 'alignment'))))
		pymol.cmd.transform_selection(tmpcutqueryname, pymolmatrix)
		pymol.cmd.align(tmpqueryname, tmpcutqueryname)
		print(tmpqueryname, tmpcutqueryname)
		pymol.cmd.delete(tmpcutqueryname)

		pymol.cmd.set_name(tmpqueryname, queryname)
		pymol.cmd.set_name(tmpsbjctname, sbjctname)

		generate_selections(obj, queryname, sbjctname)
		manage_representations(queryname, sbjctname)
		

def summary_stats(superobj, name='alignment'):
	obj = superobj['alignment']
	out = ''
	out += '\n' * 4
	out += '=' * 80 + '\n'
	out += '\n' * 2
	out += 'Summary for {}:\n'.format(name)
	out += '\tRMSD: {:0.2f} \Aa\n'.format(obj['rmsd'])
	out += '\tTM-score: {:0.2f}\n'.format(obj['quality'])
	out += '\tTransmembrane coverage: {:0.1%} of query, {:0.1%} of subject\n'.format(obj['qtmcov'], obj['stmcov'])
	out += '\tFull coverage: {:0.1%} of query, {:0.1%} of subject\n'.format(obj['qfullcov'], obj['sfullcov'])
	if 'identity' in obj:
		out += '\tIdentity in aligned region: {:0.1%}\n'.format(obj['identity'])
	return out

def rlencode(arr):
	out = []
	last = None
	for x in arr:
		if not out: out.append([x, 1])
		elif x is None and last is None: out[-1][1] += 1
		elif x is None: out.append([None, 1])
		elif last is None: out.append([x, 1])
		else: 
			if out[-1][0] + out[-1][1] == x: out[-1][1] += 1
			else: out.append([x, 1])
		last = x
	return out

def rldecode(arr):
	out = []
	for span in arr:
		if span[0] is None: out.extend([None] * span[-1])
		else: 
			try: out.extend(range(span[0], span[-1]+span[0]))
			except TypeError: out.extend([span[0]] * span[-1])
	return out

def get_resilist(pdb):
	resilist = []
	for l in pdb.split('\n'):
		if l.startswith('ATOM') and l[13:15] == 'CA': resilist.append(int(l[22:26]))
	return resilist

def compute_ranges(present, aligned, pdb):

	#actual resi
	resilist = get_resilist(pdb)
	#aligned resi
	alignmask = rldecode(aligned)

	firstpresent = present[0][0]

	print(truealigned)
	return []
	trueresi = []
	for pseudoresi in alignmask:
		if pseudoresi is None: continue
		try: trueresi.append(resilist[pseudoresi])
		#???
		except IndexError: pass
	print(trueresi)

	return trueresi
	
	

def generate_selections(superobj, qname, sname):
	obj = superobj['alignment']
	qchain = qname[5]
	schain = sname[5]


	pymol.cmd.select('qaligned', '{} and c. {} and i. {}'.format(qname, qchain, rle_to_pymol_ranges(obj['qaligned'])))
	print(obj['qaligned'])
	pymol.cmd.select('qpresent', '{} and c. {} and i. {}'.format(qname, qchain, rle_to_pymol_ranges(obj['qpresent'])))
	pymol.cmd.select('qchain', '{} and c. {}'.format(qname, qchain))

	pymol.cmd.select('saligned', '{} and c. {} and i. {}'.format(sname, schain, rle_to_pymol_ranges(obj['saligned'])))
	print(obj['saligned'])
	pymol.cmd.select('spresent', '{} and c. {} and i. {}'.format(sname, schain, rle_to_pymol_ranges(obj['spresent'])))
	pymol.cmd.select('schain', '{} and c. {}'.format(sname, schain))
	pymol.cmd.deselect()

def manage_representations(qname, sname):
	pymol.cmd.color('smudge', '{} and e. C'.format(qname))
	pymol.cmd.color('deepteal', '{} and e. C'.format(sname))
	pymol.cmd.color('splitpea', 'qchain and e. C')
	pymol.cmd.color('lightteal', 'schain and e. C')
	pymol.cmd.color('palegreen', 'qpresent and e. C')
	pymol.cmd.color('palecyan', 'spresent and e. C')
	pymol.cmd.color('tv_green', 'qaligned and e. C')
	pymol.cmd.color('cyan', 'saligned and e. C')

	pymol.cmd.hide('everything')
	pymol.cmd.show('nonbonded')
	pymol.cmd.show('lines', 'het and not r. UNK')
	pymol.cmd.show('cartoon', '*present')
	pymol.cmd.zoom('*aligned and not het')
	pymol.cmd.center('*aligned and not het')

def rle_to_pymol_ranges(arr):
	return to_pymol_ranges([[span[0], span[-1]+span[0]-1] for span in arr if span[0] is not None])

def to_pymol_ranges(arr):
	items = []
	for span in arr:
		if span[0] == span[1]: items.append(str(span[0]))
		else: items.append('{}-{}'.format(span[0], span[-1]))
	return '+'.join(items)
	
def main(args):
	if args.server: main_server(args)
	elif args.client: main_client(args)
	else: main_local(args)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	remoteparser = parser.add_mutually_exclusive_group(required=True)
	remoteparser.add_argument('--server', metavar=('DEUTDIR', 'ALIGNMENT'), nargs=2, help='Run in headless pseudoserver mode. Stdout should be accessible to a locally run instance in --client mode.')
	remoteparser.add_argument('--local', metavar=('DEUTDIR', 'ALIGNMENT'), nargs=2, help='Run in headless pseudoserver mode. Stdout should be accessible to a locally run instance in --client mode.')
	remoteparser.add_argument('--client', action='store_true', help='Run in pseudoclient mode. Must have access to stdout from quomodo running in --server mode.')

	parser.add_argument('--delay', type=float, default=2, help='How long to wait between launching Pymol and starting Pymol operations in seconds. If alignments don\'t render reliably, increasing the delay by a few seconds may help.')
	parser.add_argument('--tmalign', default='TMalign', help='TMalign path in case it doesn\'t show up on $PATH. Useful for ssh one-liners. Default: TMalign')

	args = parser.parse_args()

	main(args)
