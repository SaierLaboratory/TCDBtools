#!/usr/bin/env python

import numpy as np
from Bio import SeqIO
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-l', nargs=2, type=int, help='Return sequences between the two bounds in length')
	parser.add_argument('-m', nargs=2, type=float, help='Proportion of median')
	parser.add_argument('-z', nargs=2, type=float, help='Return sequences between the two Z-scores in length')

	parser.add_argument('infile', nargs='*')

	args = parser.parse_args()

	lengths = []
	for fn in args.infile:
		for record in SeqIO.parse(fn, 'fasta'):
			lengths.append(len(record.seq))

	if args.z:
		m = np.mean(lengths)
		s = np.std(lengths)
		bounds = [m + s * z for z in args.z]
	elif args.m:
		median = np.median(lengths)
		bounds = [median * frac for frac in args.m]
	elif args.l:
		bounds = args.l
	else: bounds = None

	i = 0
	for fn in args.infile:
		for record in SeqIO.parse(fn, 'fasta'):
			if bounds[0] <= lengths[i] <= bounds[1]: 
				print('>{}\n{}'.format(record.name, record.seq))
			i += 1


