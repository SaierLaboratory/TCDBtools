#!/usr/bin/env python

from __future__ import print_function, division
import sys
import os
import json
import argparse
import shutil
import tempfile
import subprocess

import threading
import time
UNSAFETY_DELAY = 2 #TODO: Have Pymol announce when it's done launching on its own instead of relying on this incredibly thread-unsafe thing

import pymol

import tmalignparser

def error(*stuff):
	print('[ERROR]:', *stuff, file=sys.stderr)
	sys.exit(1)


def to_pymol_ranges(spans):
	out = ''
	for span in spans:
		if span[0] == span[1]:  out += '{}+'.format(span[0])
		else: out += '{}-{}+'.format(span[0], span[-1])
	return out
		
def compress_list(l):
	lasti = None
	spans = []
	for i in l:
		if lasti is None: spans.append([i, i])
		elif i == (lasti + 1): spans[-1][1] = i
		elif i == lasti: pass
		else: spans.append([i, i])
		lasti = i

	return spans

def is_flat(l):
	if not l: return True
	elif isinstance(l[0], list): return False
	else: return True

def rle_decode(l):
	newl = []
	for x, n in l:
		newl.extend([x] * n)
	return newl

def extract_selectors(name, obj):
	query, qchain, qtms, vs, subject, schain, stms = name.split('_')
	seldict = {}


	seldict['qchain'] = '{} and c. {}'.format(query, qchain)
	seldict['schain'] = '{} and c. {}'.format(subject, schain)
	seldict['qmasked'] = '{} and c. {} and i. {}-{}'.format(
		obj['query'], 
		qchain,
		obj['qpresent'][0][0],
		obj['qpresent'][0][1]
	)
	seldict['smasked'] = '{} and c. {} and i. {}-{}'.format(
		obj['subject'], 
		schain,
		obj['spresent'][0][0],
		obj['spresent'][0][1]
	)
	alnq = []
	for span in obj['qaligned']:
		if span[0] is None: alnq += [None] * span[1]
		else: alnq += list(range(span[0], span[0]+span[1]+1))
	alns = []
	for span in obj['saligned']:
		if span[0] is None: alns += [None] * span[1]
		else: alns += list(range(span[0], span[0]+span[1]+1))
	alnd = []
	if is_flat(obj['distances']): alnd = obj['distances']
	else: alnd = rle_decode(obj['distances'])


	seldict['qpresent'] = '{} and c. {} and i. '.format(obj['query'], qchain)
	seldict['spresent'] = '{} and c. {} and i. '.format(obj['subject'], schain)

	#this will get ignored anyway
	seldict['qaligned'] = '{} and c. {} and i. '.format(obj['query'], qchain)
	seldict['saligned'] = '{} and c. {} and i. '.format(obj['subject'], schain)

	qp = obj['qpresent'][0][1] - obj['qpresent'][0][0]
	sp = obj['spresent'][0][1] - obj['spresent'][0][0]
	qa = 0
	sa = 0

	qpresent = []
	spresent = []

	for q, d, s in zip(alnq, alnd, alns):
		if q is not None: qpresent.append(q)
		if s is not None: spresent.append(s)
		if (q is not None) and (s is not None) and (d is not None):
			qpresent.append(q)
			spresent.append(s)
			qa += 1
			sa += 1

	seldict['qpresent'] += to_pymol_ranges(compress_list(qpresent))
	seldict['spresent'] += to_pymol_ranges(compress_list(spresent))

	print('='*80)
	print('{}:'.format(name))
	print('\tRMSD: {:0.2f}'.format(obj['rmsd']))
	print('\tTM-score: {:0.2f}'.format(obj['quality']))
	if qa/qp < sa/sp:
		print('\tMinimum coverage: {:0.2%} of query'.format(qa/qp))
		print('\tMaximum coverage: {:0.2%} of subject'.format(sa/sp))
	else:
		print('\tMinimum coverage: {:0.2%} of subject'.format(sa/sp))
		print('\tMaximum coverage: {:0.2%} of query'.format(qa/qp))
	seldict['qpresent'] = seldict['qpresent'][:-1]
	seldict['spresent'] = seldict['spresent'][:-1]
	seldict['qaligned'] = seldict['qaligned'][:-1]
	seldict['saligned'] = seldict['saligned'][:-1]
	return seldict

def run_tmap(obj, d2dir):
	qfn = '{}/cut_pdbs/{}'.format(d2dir, os.path.basename(obj['queryfn']))
	sfn = '{}/cut_pdbs/{}'.format(d2dir, os.path.basename(obj['subjectfn']))

	tf = tempfile.NamedTemporaryFile()
	out = subprocess.check_output(['TMalign', qfn, sfn])
	tf.write(out)
	tf.flush()

	selectors = tmalignparser.main([tf.name, qfn, sfn])
	return selectors


def main(name, obj, d2dir, wd=None, outfile=None):


	if wd is None:
		print('Error: No output directory specified')
		exit(1)

	pymol.cmd.set('fetch_path', wd)
	#pymol.cmd.fetch(obj['query'])
	pymol.cmd.load('{}/pdbs/{}.pdb'.format(d2dir, obj['query']))
	#pymol.cmd.fetch(obj['subject'])
	pymol.cmd.load('{}/pdbs/{}.pdb'.format(d2dir, obj['subject']))

	selectors = extract_selectors(name, obj)

	qaligned, saligned = run_tmap(obj, d2dir)
	selectors['qaligned'] = qaligned
	selectors['saligned'] = saligned

	pymol.cmd.hide('lines')
	pymol.cmd.hide('spheres', 'r. DUM')
	pymol.cmd.show('nonbonded', 'r. DUM')
	print()
	print('Selectors:')
	print('----------')
	for k in sorted(selectors):
		if k != 'matrix':
			pymol.cmd.select(k, selectors[k])
			
			if k.endswith('masked'): 
				pymol.cmd.show('cartoon', k)
			print(k + ':', selectors[k])


	truematrix = []
	for row in obj['matrix']: truematrix += row
	truematrix += [0., 0., 0., 1.]
	pymol.cmd.transform_selection(obj['query'], truematrix)
	pymol.cmd.rotate('x', 90)
	pymol.cmd.center('*present')
	pymol.cmd.zoom('*present')
	pymol.cmd.show('cartoon', '*present')

	#pymol.cmd.color('palegreen', obj['query'])
	#pymol.cmd.color('palecyan', obj['subject'])
	pymol.cmd.color('smudge', obj['query'])
	pymol.cmd.color('deepteal', obj['subject'])
	pymol.cmd.color('splitpea', 'qchain')
	pymol.cmd.color('lightteal', 'schain')
	pymol.cmd.color('palegreen', 'qmasked')
	pymol.cmd.color('palecyan', 'smasked')
	pymol.cmd.color('tv_green', 'qaligned')
	pymol.cmd.color('cyan', 'saligned')
	#pymol.cmd.remove('r. DUM')

	pymol.cmd.deselect()
	pymol.cmd.rotate('z', 360)


	if outfile: 
		pymol.cmd.save(outfile)

		#maybe work in the silent option instead of autoquitting
		pymol.cmd.quit()

class DoStuffThread(threading.Thread):
	def run(self):
		time.sleep(UNSAFETY_DELAY)
		main(**self._kwargs)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('d1dir', help='Deuterocol1 directory (or deuterocol files directory in 2-arg mode)')
	parser.add_argument('d2dir', help='Deuterocol2 directory (or alignment name in 2-arg mode)')

	parser.add_argument('name', nargs='?', help="name of alignment to load. Be warned that this selects only the first occurrence in all the files in the order given. If this argument is left out, `d1dir' is interpreted as an all-containing deuterocol.py-generated directory, and `d2dir' is interpreted as `name'")

	parser.add_argument('-o', default=None, help='Save the resulting session')

	args = parser.parse_args()

	if not args.name:
		d1dir = args.d1dir + '/deuterocol1'
		d2dir = args.d1dir
		name = args.d2dir
	else: 
		d1dir = args.d1dir
		d2dir = args.d2dir
		name = args.name
		

	query, qchain, qhel, vs, subject, schain, shel = name.split('_')
	possfams = set()

	with open(d1dir + '/pdblist.json') as f: pdbmap = json.loads(f.read())
	for fam in pdbmap:
		if (query + '_' + qchain) in pdbmap[fam]: possfams.add(fam)
		if (subject + '_' + schain) in pdbmap[fam]: possfams.add(fam)

	possfams = tuple(possfams)
	subdirlist = []
	for subdir in os.listdir(d2dir):
		count = 0
		for fam in possfams: 
			if fam in subdir: count += 1
		if count >= 2: subdirlist.append(subdir)
	if not subdirlist: error('Could not find directories for comparisons between families {}'.format(possfams))

	found = False
	jstr = ''
	n = 0

	for subdir in subdirlist:
		fn = d2dir + '/' + subdir + '/tmalignments/sp_all.tsv'
		if not os.path.isfile(fn): error('Could not find TMalign records for', subdir)

	#for fn in args.infile:
		with open(fn) as f:
			for l in f:
				n += 1
				if not n % 10000: print('Checked', n, 'lines')
				if l.startswith(name[:2]):
					if l.startswith(name):
						print('Found', name, 'in', fn)
						jstr = l[l.find('\t')+1:]
						found = True
				if found: break
		if found: break

	if not jstr: 
		print('Error: Could not find {} in {}'.format(name, infile))

	obj = json.loads(jstr)
	wd = tempfile.mkdtemp()

	#implement headless
	sys.argv = sys.argv[:1]
	if args.o: pymol.pymol_argv = ['pymol','-qc']

	try: 
		dsthread = DoStuffThread(kwargs={'name':name, 'obj':obj, 'wd':wd, 'd2dir':d2dir, 'outfile':args.o})
		dsthread.start()
		pymol.launch()
	finally: shutil.rmtree(wd)
