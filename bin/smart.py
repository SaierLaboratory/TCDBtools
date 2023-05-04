#!/usr/bin/env python3
#-*- encoding: utf-8 -*-
#stochastic multiple alignment refinement tool
from __future__ import print_function, division

from Bio import SeqIO, AlignIO
import numpy as np

import subprocess
import argparse
import tempfile
import sys
import time
import io
import os

VERBOSITY = 0
GREETING = [1]
def debug(*things):
	if VERBOSITY < 3: return
	matters = [str(x) for x in things]
	print('Foreman: "Have you considered buying new equipment, sir?', ' '.join(matters) + '"',file=sys.stderr)
def info(*things): 
	if VERBOSITY < 2: return
	matters = [str(x) for x in things]
	if GREETING: 
		t = time.localtime()
		if t.tm_hour <= 4: greetstr = 'Good evening, sir!'
		elif t.tm_hour < 12: greetstr = 'Good morning, sir!'
		elif t.tm_hour < 18: greetstr = 'Good afternoon, sir!'
		elif t.tm_hour < 24: greetstr = 'Good evening, sir!'
		else: greetstr = 'Top of the day to you, sir!'
		print('Foreman: "' + greetstr, ' '.join(matters) + '"',file=sys.stderr)
		GREETING.pop()
	else: 
		r = np.random.randint(1, 5)
		if r <= 1: introstr = '''There's not much to report, sir.'''
		elif r <= 2: introstr = '''Everything's in order, sir.'''
		elif r <= 3: introstr = '''Work continues, sir.'''
		elif r <= 4: introstr = '''Progress continues.'''
		elif r <= 5: introstr = '''Ah, the sounds of work in the foundry!'''
		print('Foreman: "' + introstr, ' '.join(matters) + '"',file=sys.stderr)
def apologize(*things): 
	if VERBOSITY < 1: return
	matters = [str(x) for x in things]
	print('Foreman: "Apologies, sir.', ' '.join(matters) + '"',file=sys.stderr)
def warn(*things): 
	if VERBOSITY < 1: return
	matters = [str(x) for x in things]
	print('Foreman: "Do be careful, sir.', ' '.join(matters) + '"',file=sys.stderr)
def error(*things): 
	matters = [str(x) for x in things]
	print('Foreman: "There\'s been an accident in the foundry!', ' '.join(matters) + '"',file=sys.stderr)
	exit(1)
def critical(*things): 
	matters = [str(x) for x in things]
	print('Foreman: "The inspector\'s coming for our books!', ' '.join(matters) + '"',file=sys.stderr)
	exit(1)

class Refinery(object):
	def __init__(self, msa):
		self.msa = msa
		self.orig_msa = msa

	def get_occupancies(self):
		columns = np.zeros(self.msa.get_alignment_length())
		for seq in self.msa:
			for i, c in enumerate(seq.seq):
				if c != '-': columns[i] += 1
		return columns / len(self.msa)
		
	def get_median_occupancy(self): return np.median(self.get_occupancies())

	def get_mean_occupancy(self): return np.mean(self.get_occupancies())

	def get_num_sequences(self): return len(self.msa)

	def get_alignment_length(self): return self.msa.get_alignment_length()

	def run(self):
		print(self.get_median_occupancy())

	def sift(self, target=0.5, write_off=0.2):
		occ = self.get_median_occupancy()
		info('The raw ore is {:0.0%} pure! We are separating out the obvious impurities as we speak.'.format(occ))

		results = {}
		fmsa = io.StringIO()
		AlignIO.write(self.msa, fmsa, 'clustal')
		fmsa.flush()
		fmsa.seek(0)

		tfmsa = tempfile.NamedTemporaryFile()
		AlignIO.write(self.msa, tfmsa.name, 'clustal')
		tfmsa.flush()

		for resoverlap in np.arange(0.5, 0.90+0.1, 0.1):
			for seqoverlap in np.arange(50, 90+10, 10):
				
				cmd = ['trimal', '-resoverlap', str(resoverlap), '-seqoverlap', str(seqoverlap), '-in', tfmsa.name]

				p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

				out, err = p.communicate()

				if not out.decode('utf-8').strip(): continue

				fout = io.StringIO()
				fout.write(out.decode('utf-8'))
				fout.flush()
				fout.seek(0)
				newmsa = AlignIO.read(fout, 'clustal')
				results[(resoverlap, seqoverlap)] = Refinery(newmsa)

		sortme = []
		for k in results:
			damage = (self.get_num_sequences() - results[k].get_num_sequences()) / self.get_num_sequences()
			occupancy = results[k].get_median_occupancy()
			sortme.append((occupancy, -damage, k, results[k]))

		bestresult = None
		for row in reversed(sorted(sortme)):
			if row[1] > write_off: continue
			else:
				bestresult = row[-1]
				info('We were able to sift the ore to a median purity of {:0.0%}, sir! We used a resoverlap of {:0.0%} and a seqoverlap of {}%, throwing out {:0.0%} as insufficiently sifted.'.format(row[0], row[2][0], row[2][1], -row[1]))
				if row[0] < target: apologize('That falls somewhat short of your target of {:0.0%}.'.format(target))
				break

		if bestresult is None: error('The ore could not be sifted!')
		return bestresult

	def enrich(self, target=0.5, write_off=0.1):
		info('We are now enriching the remnants of the ore to remove unwanted rock.')

		results = {}
		fmsa = io.StringIO()
		AlignIO.write(self.msa, fmsa, 'clustal')
		fmsa.flush()
		fmsa.seek(0)

		tfmsa = tempfile.NamedTemporaryFile()
		AlignIO.write(self.msa, tfmsa.name, 'clustal')
		tfmsa.flush()

		for gapthreshold in np.arange(0.05, 0.90+0.05, 0.05):
			
			cmd = ['trimal', '-terminalonly', '-gt', str(gapthreshold), '-in', tfmsa.name]

			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

			out, err = p.communicate()

			if not out.decode('utf-8').strip(): continue

			fout = io.StringIO()
			fout.write(out.decode('utf-8'))
			fout.flush()
			fout.seek(0)
			newmsa = AlignIO.read(fout, 'clustal')
			results[gapthreshold] = Refinery(newmsa)

		sortme = []
		for k in results:
			damage = (self.get_alignment_length() - results[k].get_alignment_length()) / self.get_alignment_length()
			occupancy = results[k].get_mean_occupancy()
			sortme.append((occupancy, -damage, k, results[k]))

		bestresult = None
		for row in reversed(sorted(sortme)):
			if -row[1] > write_off: continue
			else:
				bestresult = row[-1]
				info('We were able to enrich the ore to a median purity of {:0.0%}, sir! We used a gapthreshold of {:0.0%}, throwing out {:0.0%} as insufficiently pure.'.format(row[0], row[2], -row[1]))
				if row[0] < target: apologize('That falls somewhat short of your target of {:0.0%}.'.format(target))
				break

		if bestresult is None: error('The ore could not be enriched without damaging more than {:0.0%}!'.format(write_off))
		return bestresult

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='?', default='/dev/stdin', help='File to read in (default: stdin)')
	parser.add_argument('-o', '--outdir', required=True, help='Where to store the refined ore')

	horizontalargs = parser.add_argument_group('horizontal refinement arguments')
	horizontalargs.add_argument('-hd', '--acceptable-horizontal-collateral-damage', type=float, default=0.2, help='Maximum fraction of sequences to eliminate during the horizontal refinement step (default: 0.2)')
	horizontalargs.add_argument('-ho', '--target-median-occupancy', type=float, default=0.5, help='Target minimum occupancy at each column')

	verticalargs = parser.add_argument_group('vertical refinement arguments')
	verticalargs.add_argument('-vd', '--acceptable-vertical-collateral-damage', type=float, default=0.5, help='Maximum net reduction in alignment length during the terminalonly vertical refinement step (default: 0.5)')
	verticalargs.add_argument('-vo', '--target-mean-occupancy', type=float, default=0.4, help='Target minimum occupancy at each column')

	args = parser.parse_args()

	if not os.path.isdir(args.outdir): os.mkdir(args.outdir)
	with open('{}/commandline'.format(args.outdir), 'w') as fh:
		fh.write('''### SEQUENCEWISE PARAMS ###
# target: {}
# max damage: {}

### COLUMNWISE PARAMS ###
# target: {}
# max damage: {}

{}'''.format(args.target_median_occupancy, args.acceptable_horizontal_collateral_damage, args.target_mean_occupancy, args.acceptable_vertical_collateral_damage, ' '.join(sys.argv)))

	VERBOSITY = 2
	info('We are loading the ores and fuel.')
	with open(args.infile) as f:
		for msa in AlignIO.parse(f, 'clustal'):
			refinery = Refinery(msa)
			AlignIO.write(refinery.msa, '{}/original.aln'.format(args.outdir), 'clustal')
			result = refinery.sift(target=args.target_median_occupancy, write_off=args.acceptable_horizontal_collateral_damage)
			refinery.msa = result.msa
			AlignIO.write(refinery.msa, '{}/step1.aln'.format(args.outdir), 'clustal')
			
			result = refinery.enrich(target=args.target_mean_occupancy, write_off=args.acceptable_vertical_collateral_damage)
			refinery.msa = result.msa
			AlignIO.write(refinery.msa, '{}/step2.aln'.format(args.outdir), 'clustal')
