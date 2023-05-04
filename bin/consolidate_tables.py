#!/usr/bin/env python

from kevtools_common.types import TCID
import argparse
import os
import json2tsv

def main(tabledir, outdir, level=3):
	if not os.path.isdir(outdir): os.mkdir(outdir)
	infhs = {}
	for bn in os.listdir(tabledir):
		if '_vs_' in bn and '.tsv' in bn:
			famvfam = bn[:bn.find('.tsv')]
			fams = famvfam.split('_vs_')
			fam1 = TCID(fams[0])[:level]
			fam2 = TCID(fams[1])[:level]

			if (fam1, fam2) not in infhs: infhs[fam1,fam2] = []
			infhs[fam1,fam2].append(open('{}/{}'.format(tabledir, bn)))
	if not infhs and os.path.isdir('{}/tables'.format(tabledir)):
		for bn in os.listdir('{}/tables'.format(tabledir)):
			if '_vs_' in bn and '.tsv' in bn:
				famvfam = bn[:bn.find('.tsv')]
				fams = famvfam.split('_vs_')
				fam1 = TCID(fams[0])[:level]
				fam2 = TCID(fams[1])[:level]

				if (fam1, fam2) not in infhs: infhs[fam1,fam2] = []
				infhs[fam1,fam2].append(open('{}/tables/{}'.format(tabledir, bn)))

	outfhs = {}
	for pair in infhs:
		outfhs[pair] = open('{}/{}_vs_{}.tsv'.format(outdir, pair[0], pair[1]), 'w')

	for fampair in infhs:
		outfile = outfhs[fampair]
		for infile in infhs[fampair]:
			outfile.write(infile.read())

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--level', type=int, default=3, help='TC-level to consolidate table results to. (default:3)')

	parser.add_argument('tabledir', help='Directory containing tables to concatenate. Files must be named [TCID]_vs_[TCID].tsv')
	parser.add_argument('outdir', help='Directory to write results to')

	args = parser.parse_args()

	main(args.tabledir, args.outdir, level=args.level)
