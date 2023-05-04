#!/usr/bin/env python
from Bio import SeqIO
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='*', default=['/dev/stdin'])

	args = parser.parse_args()

	for fn in args.infile:
		with open(fn) as fh:
			for record in SeqIO.parse(fh, 'fasta'):
				print('{}\t{}'.format(record.name, len(record.seq)))
