#!/usr/bin/env python

import os
import sys
import argparse
from Bio import SeqIO

def warn(*things): print('[WARNING]', *things, file=sys.stderr)
def error(*things): 
	print('[ERROR]', *things, file=sys.stderr)
	exit(1)

def parse_cdhit(fh):
	clusters = {}
	clname = None
	for l in fh:
		if l.startswith('>'):
			clname = l[1:-1]
			clusters[clname] = []

		else:
			acc10 = l[l.find('>')+1:l.find('...')]
			clusters[clname].append(acc10)

	return clusters

def rev_dict(d): 
	newdict = {}
	for k in d: 
		for v in d[k]: newdict[v] = k
	return newdict

def main(clstrfn, fastafn, outdir, ignore_missing=False):

	with open(clstrfn) as fh: clusters = rev_dict(parse_cdhit(fh))

	by_cluster = {}
	with open(fastafn) as fh:
		for record in SeqIO.parse(fh, 'fasta'):
			if record.id not in clusters: 
				if ignore_missing:
					warn('Somehow, {} is missing!'.format(record.id))
					continue
				else:
					error('Somehow, {} is missing!'.format(record.id))
					
			if clusters[record.id] not in by_cluster: by_cluster[clusters[record.id]] = []
			by_cluster[clusters[record.id]].append(record)

	if not os.path.isdir(outdir): os.mkdir(outdir)
	for clname in by_cluster:
		clean_clname = clname.replace(' ', '_').replace('/', '_')
		SeqIO.write(by_cluster[clname], '{}/{}.faa'.format(outdir, clean_clname), 'fasta')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('clstrfn', help='.clstr file CD-HIT outputs')
	parser.add_argument('fastafn', help='Complete FASTA with all sequences (input to CD-HIT')

	parser.add_argument('-o', required=True, help='Where to write the resulting FASTAs')

	parser.add_argument('--ignore-missing', action='store_true', help='Don\'t crash on missing sequencse')
	args = parser.parse_args()

	main(args.clstrfn, args.fastafn, args.o, ignore_missing=args.ignore_missing)
