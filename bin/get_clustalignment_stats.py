#!/usr/bin/env python

from Bio import AlignIO
import smart
import argparse
import numpy as np

def dict_entropy(d):
	total = sum(d.values())
	shannon = 0
	for k in d:
		p = d[k]/total
		shannon -= p * np.log2(p)
	return shannon

def get_msa_entropy(refinery):
	entropies = []
	pos = {}
	for seq in refinery.msa:
		for i, r in enumerate(seq):
			if i not in pos: pos[i] = {}
			if r not in pos[i]: pos[i][r] = 0
			pos[i][r] += 1
	for i in sorted(pos):
		entropies.append(dict_entropy(pos[i]))

	return entropies

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', default=['/dev/stdin'], nargs='*', help='MSAs to read in')

	args = parser.parse_args()

	for fn in args.infile:
		for msa in AlignIO.parse(fn, 'clustal'):
			refinery = smart.Refinery(msa)
			print(fn)
			print('\t{}x{} alignment'.format(refinery.get_num_sequences(), refinery.get_alignment_length()))

			occ = refinery.get_occupancies()
			print('\tMedian column occupancy: {:0.1%}'.format(np.median(occ)))
			print('\tMean column occupancy: {:0.1%} +/- {:0.1%}'.format(np.mean(occ), np.std(occ)))
			ent = get_msa_entropy(refinery)
			print('\tMedian column heterogeneity: {:0.1%}'.format(np.median(ent)/np.log2(20)))
			print('\tMean column heterogeneity: {:0.1%} +/- {:0.1%}'.format(np.mean(ent)/np.log2(20), np.std(ent)/np.log2(20)))

			print()

