#!/usr/bin/env python

import argparse
import kevtools_common.types
import numpy as np

def import_alignments(fh, expect_cutoff=1.0, tclevel=3):
	if isinstance(fh, str): fh = open(fh)

	alignments = {}

	for l in fh:
		if not l.strip(): continue
		elif l.strip().startswith('#'): continue
		else:
			sl = l.replace('\n', '').split('\t')

			qtcid = kevtools_common.types.TCID(sl[0])
			stcid = kevtools_common.types.TCID(sl[1])
			pident = float(sl[2])
			qspan = (int(sl[6]), int(sl[7]))
			sspan = (int(sl[7]), int(sl[8]))
			expect = float(sl[10])
			bits = float(sl[11])
			qlen = int(sl[12])
			slen = int(sl[13])

			if expect > expect_cutoff: continue

			#defined as families with >200 sequences
			#bigfams = ('3.A.1', '2.A.1', '2.A.7', '3.D.1', '3.A.7', '1.A.1', '3.A.2', '3.A.3', '2.A.6')
			#mod = 1 if str(qtcid[:3]) in bigfams else 0
			mod = 0

			if qtcid[:tclevel+mod] not in alignments: alignments[qtcid[:tclevel+mod]] = []
			#if qtcid[:3] not in alignments: alignments[qtcid[:3]] = {}
			#if stcid[:3] not in alignments: alignments[stcid[:3]] = {}

			#if stcid[:3] not in alignments[qtcid[:3]]: alignments[qtcid[:3]][stcid[:3]] = []
			#if qtcid[:3] not in alignments[stcid[:3]]: alignments[stcid[:3]][qtcid[:3]] = []

			#alignments[qtcid[:3]][stcid[:3]].append((expect, qtcid, stcid))
			#alignments[stcid[:3]][qtcid[:3]].append((expect, stcid, qtcid))

			obj = {'query':qtcid, 'subject':stcid,
				'pident':pident,
				'expect':expect, 'bits':bits,
				'qspan':qspan, 'sspan':sspan,
				'qlen':qlen, 'slen':slen}

			#alignments[qtcid[:3]].append((bits, expect, qtcid, stcid))
			alignments[qtcid[:tclevel+mod]].append(obj)

	return alignments

def get_center(alnlist, param='bits'):
	totbits = {}
	for aln in alnlist:
		if aln['query'] != aln['subject']: 
			if aln['query'] not in totbits:
				totbits[aln['query']] = 0
			totbits[aln['query']] += aln[param]
	if not totbits: return None

	keys, values = totbits.keys(), totbits.values()
	valarr = sorted(values)
	#print(max(valarr), np.mean(valarr), '+/-', np.std(valarr))

	maxbits = max(values)
	for k, v in zip(keys, values):
		if v == maxbits: return k

class ScoreFunction(object):
	@staticmethod
	def get_method(name):
		self = ScoreFunction
		if name.startswith('fixed'): return self.fixed
		elif name.startswith('bits'): return self.bits
		elif name.startswith('sqrt_bits'): return self.sqrt_bits
		elif name.startswith('avg_length'): return self.avg_length
		elif name.startswith('sqrt_avg_length'): return self.sqrt_avg_length
		elif name.startswith('pident'): return self.pident
		else: raise AttributeError('No matching method for `{}\''.format(name))

	@staticmethod
	def fixed(obj): return 1, 1

	#FIXME: two returned, not one
	@staticmethod
	def bits(obj): return obj['bits']

	@staticmethod
	def sqrt_bits(obj): return obj['bits']**.5

	@staticmethod
	def avg_length(obj): return np.mean([obj[k][1] - obj[k][0] + 1 for k in ('qspan', 'sspan')])

	@staticmethod
	def sqrt_avg_length(obj): return np.mean([obj[k][1] - obj[k][0] + 1 for k in ('qspan', 'sspan')])**.5

	@staticmethod
	def pident(obj): return obj['pident'] * 0.01

def calculate_goodspans(alnlist, scorefunction=ScoreFunction.fixed):
	seqdict = {}
	for aln in alnlist:
		if aln['query'] not in seqdict: seqdict[aln['query']] = np.zeros(aln['qlen'])
		if aln['subject'] not in seqdict: seqdict[aln['subject']] = np.zeros(aln['slen'])

		value = scorefunction(aln)
		for i in range(aln['qspan'][0]-1, aln['qspan'][1]):
			seqdict[aln['query']][i] += value
		for i in range(aln['sspan'][0]-1, aln['sspan'][1]):
			seqdict[aln['subject']][i] += value

	return seqdict

def prettyprint(scorevec):
	#norm = goodspans[tcid] / max(goodspans[tcid])
	norm = scorevec

	out = ''
	for c in norm:
		if c <= 0.3: out += '.'
		elif c <= 0.7: out += ':'
		else: out += '|'
	return out

def numerate(seq, split=80, add_ladder=True):
	out = ''

	if add_ladder: 
		ladder = ''
		for i in range(10, split+1, 10): ladder += str(i).rjust(10)
		out += '#' + ladder[1:]
		out += '\n'

	for i in range(0, len(seq), split):
		out += seq[i:i+split]
		out += '\n'
	return out.rstrip()

def parse_intrafam_output(fh):
	# TODO: FIXME: IMPLEMENT ME
	with open(fh) as f:
		tcid = None
		fam = None
		results = {}
		for l in f:
			if l.startswith('>>'): 
				tcid = kevtools_common.types.TCID(l[1:].strip())
				results[fam].something
			elif l.startswith('>'):
				fam = 0
	pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='?', default='/dev/stdin')
	parser.add_argument('--evalue', default=1e-3, type=float, help='Maximum expect value to allow between family "members"')
	args = parser.parse_args()

	with open(args.infile) as fh:
		alignments = import_alignments(fh, expect_cutoff=args.evalue)
		for fam in alignments:

			scorefunction = ScoreFunction.fixed

			#print('{}\t{}'.format(fam, get_center(alignments[fam])))
			goodspans = calculate_goodspans(alignments[fam], scorefunction=scorefunction)
			normalized = {}
			#totscore = 0
			#for aln in alignments[fam]: totscore += scorefunction(aln)
			for k in goodspans: 
				normalized[k] = goodspans[k] / max(goodspans[k])

			print('>{} ({} members)'.format(fam, len(goodspans)))
			for tcid in goodspans:
				
				print('>>{}'.format(tcid))
				print('# Raw max: {}'.format(max(goodspans[tcid])))
				print('{}\n'.format(numerate(prettyprint(normalized[tcid]))))
