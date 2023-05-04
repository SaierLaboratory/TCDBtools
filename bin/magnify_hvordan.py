#!/usr/bin/env python

from __future__ import print_function, division

import argparse
import os
import subprocess
import re
import shutil

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec

import quod

#python3 note: capitalization will change slightly, but py2to3 should take care of it
#import ConfigParser

class Plot(quod.Entity):
	def __init__(self):
		self.entities = []

class Protocol2(object):
	def __init__(self, p2dir):
		self.p2dir = p2dir

	def get_hit(self, pair=None):
		try: iter(pair)
		except TypeError: raise TypeError('pair must be iterable')
		def parse_for(fh, accs):
			for l in fh:
				if not l.strip(): continue
				elif l.startswith('#'): continue
				sl = l.strip().split('\t')
				subject, target, ssearch, gsat, slen, tlen, sseq, tseq, tmo = sl[:9]
				if subject in accs and target in accs: 
					return {'subject':subject, 'target':target, 'ssearch':ssearch, 'gsat':float(gsat), 'slen':int(slen), 'tlen':int(tlen), 'sseq':sseq, 'tseq':tseq, 'tmo':float(tmo)}
			return None

		#TODO: unified proto2 traversal class
		path = os.getcwd()
		path += '/{}'.format(self.p2dir)
		found = False
		for subdir in os.listdir(path):
			path += '/{}'.format(subdir)
			for bn in os.listdir(path):
				if bn == 'report.tbl':
					#found = True
					fn = path + '/{}'.format(bn)
					with open(fn) as f:
						found = parse_for(f, pair)
				if found: 
					break
			if found: 
				found['fams'] = subdir.split('_vs_')
				break
		if not found: raise IOError('Could not find {} and {} in report.tbl'.format(*pair))

		return found
#class Protocol1(object):
#	def __init__(self, p1dir):
#		self.p1dir = p1dir
#
#	def get_hit(self, acc)

#TODO: allow opening psiblast.tbl directly
def get_p1hit(p1dir, accs):
	def parse_for(fh, acclist):
		acclistcopy = acclist[:]
		found = {}
		for l in fh:
			if not acclistcopy: break
			if not l.strip(): continue
			elif l.startswith('#'): continue
			sl = l.strip().split('\t')
			query, sacc, subject, score, evalue, pident, qstart, qend, qlen, sstart, send, slen, ssciname, qseq, sseq = sl[:15]
			for i, acc in enumerate(acclistcopy):
				if acc in l:
					acclistcopy.pop(i)
					found[acc] = {}
					found[acc]['query'] = query
					found[acc]['sacc'] = sacc 
					found[acc]['subject'] = subject
					found[acc]['score'] = float(score)
					found[acc]['evalue'] = float(evalue)
					found[acc]['pident'] = float(pident)
					found[acc]['qstart'] = int(qstart)
					found[acc]['qend'] = int(qend)
					found[acc]['qlen'] = int(qlen)
					found[acc]['sstart'] = int(sstart)
					found[acc]['send'] = int(send)
					found[acc]['slen'] = int(slen)
					found[acc]['ssciname'] = ssciname
					found[acc]['qseq'] = qseq
					found[acc]['sseq'] = sseq
				if acc in found: continue
		return found
									
	found = {}
	with open(p1dir + '/psiblast.tbl') as f:
		found.update(parse_for(f, accs))
	return found

def fasta_length(fastastr):
	if not fastastr.startswith('>'): return len(re.sub('\s', '', fastastr))
	else: return len(re.sub('\s', '', fastastr[fastastr.find('\n'):]))

def main(fxdir, p2dir, acc1, acc2, outdir='magnify_plots'):
	p2 = Protocol2(p2dir=p2dir)
	accs = [acc1, acc2]
	p2hit = p2.get_hit(accs)

	p1hits = {}
	for fam in p2hit['fams']:
		path = fxdir + '/' + fam
		p1hits.update(get_p1hit(path, [acc1, acc2]))
		#p1 = Protocol1(p1dir=fxdir)
	if not os.path.isdir(outdir): os.mkdir(outdir)
	if not os.path.isdir(outdir + '/sequences'): os.mkdir(outdir + '/sequences')
	tcaccs = []
	for acc in accs:
		if not os.path.isfile(outdir + '/sequences/' + acc + '.fa'):
			seq = subprocess.check_output(['blastdbcmd', '-db', 'nr', '-target_only', '-entry', acc])
			with open(outdir + '/sequences/' + acc + '.fa', 'w') as f: f.write(seq)

		if not os.path.isfile(outdir + '/sequences/' + p1hits[acc]['query'] + '.fa'):
			seq = subprocess.check_output(['blastdbcmd', '-db', 'tcdb', '-target_only', '-entry', p1hits[acc]['query']])
			with open(outdir + '/sequences/' + p1hits[acc]['query'] + '.fa', 'w') as f: f.write(seq)

		if not os.path.isfile(outdir + '/sequences/' + 'p1-' + p1hits[acc]['query'] + '.fa'):
			seq = '>' + p1hits[acc]['query'] + '\n' + p1hits[acc]['qseq']
			with open(outdir + '/sequences/' + 'p1-' + p1hits[acc]['query'] + '.fa', 'w') as f: f.write(seq)

		if not os.path.isfile(outdir + '/sequences/' + 'p1-' + p1hits[acc]['sacc'] + '.fa'):
			seq = '>' + p1hits[acc]['sacc'] + '\n' + p1hits[acc]['sseq']
			with open(outdir + '/sequences/' + 'p1-' + p1hits[acc]['sacc'] + '.fa', 'w') as f: f.write(seq)

	if not os.path.isfile(outdir + '/sequences/' + 'p2-' + p2hit['subject'] + '.fa'):
		seq = '>' + p2hit['subject'] + '\n' + p2hit['sseq']
		with open(outdir + '/sequences/' + 'p2-' + p2hit['subject'] + '.fa', 'w') as f: f.write(seq)

	if not os.path.isfile(outdir + '/sequences/' + 'p2-' + p2hit['target'] + '.fa'):
		seq = '>' + p2hit['target'] + '\n' + p2hit['tseq']
		with open(outdir + '/sequences/' + 'p2-' + p2hit['target'] + '.fa', 'w') as f: f.write(seq)

	#FIXME: configurability
	ab_mult = 1.
	cd_mult = 1.
	ef_mult = 1.
	g_mult = 0.5


	fig = Figure()
	canvas = FigureCanvas(fig)

	gs_ab = gridspec.GridSpec(2, 1)
	gs_ab.update(top=0.02, bottom=0.23)
	ax_a  = fig.add_subplot(121)
	ax_b  = fig.add_subplot(122)
	print(dir(gs_ab))

	#gs_cd.update(top=0.27, bottom=0.48)
	#ax_c  = fig.add_subplot(121)
	#ax_d  = fig.add_subplot(122)
	fig.savefig('test.png')
	#gs_cd = gridspec.GridSpec(2, 1)
	#gs_ef = gridspec.GridSpec(2, 1)
	#gs_g = gridspec.GridSpec(2, 1)

	#fig.savefig('test.png')

	if not os.path.isdir(outdir + '/config'): os.mkdir(outdir + '/config')
	
	
	#print('p1', p1hits)
	#print('p2', p2hit)

	pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-n', action='store_true', help='just generate a configureation file and exit')

	parser.add_argument('-o', '--outdir', default='magnify_plots', help='output directory')

	parser.add_argument('fxdir', help='famXpander directory')
	parser.add_argument('p2dir', help='Protocol2 directory')
	#IMPLEMENT direct plotting from HTML
	#parser.add_argument('--mode'

	parser.add_argument('acc1', help='B sequence')
	parser.add_argument('acc2', help='C sequence')

	args = parser.parse_args()

	#main(args.fxdir, args.p2dir, [[args.acc1, args.acc2]])
	

	main(args.fxdir, args.p2dir, args.acc1, args.acc2, outdir=args.outdir)
