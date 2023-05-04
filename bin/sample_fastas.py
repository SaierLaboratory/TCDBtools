#!/usr/bin/env python

from Bio import SeqIO
import random
import argparse
import os
from kevtools_common import messages

#print('{:016x}'.format(seed))

def main(infnlist, samplesize, samples=1, outdir=None, initseed=0, replacement=False, merge=False):
	seed = initseed

	if merge:
		pass
	else:
		if outdir is None: outdir = 'samples_{:016x}'.format(abs(hash(tuple(infnlist))))
		if not os.path.isdir(outdir): os.mkdir(outdir)

	seqcount = 0
	for fn in infnlist:
		if fn == '/dev/stdin': messages.info('Now reading from stdin')
		with open(fn) as f:
			for record in SeqIO.parse(f, 'fasta'):
				seqcount += 1

	if not replacement:
		pickme = list(range(seqcount))
		random.shuffle(pickme)

	if merge: open(outdir, 'w')

	for sample in range(samples):
		random.seed(seed)
		#outfile = '{}{:016x}.faa'.format(prefix, seed)

		if replacement:
			pickme = list(range(seqcount))
			random.shuffle(pickme)
			pickme2 = set(pickme[:samplesize])
		else: 
			pickme2 = set(pickme[sample*samplesize:(sample+1)*samplesize])
			if not pickme2: break

		if pickme2: lasti = max(pickme)
		else: lasti = None

		if merge:
			i = 0
			with open(outdir, 'a') as outf:
				for fn in infnlist:
					with open(fn) as f:
						for record in SeqIO.parse(f, 'fasta'):
							if i in pickme2: outf.write('>{}\n{}\n'.format(record.name, record.seq))
							elif i > lasti: break
							i += 1
						
		else:
			i = 0
			with open('{}/sample{}.faa'.format(outdir, sample), 'w') as outf:
				for fn in infnlist:
					with open(fn) as f:
						for record in SeqIO.parse(f, 'fasta'):
							if i in pickme2: outf.write('>{}\n{}\n'.format(record.name, record.seq))
							elif i > lasti: break
							i += 1

		if replacement: seed = random.randint(0, 0xffffffffffffffff)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='*', default=['/dev/stdin'], help='FASTAs to sample from (default: stdin)')
	parser.add_argument('-c', metavar='SAMPLESIZE', type=int, required=True, help='How many sequences to place in each sample')
	parser.add_argument('-m', '--merge', action='store_true', help='Concatenate sequences and write them to the \033[1mfile\033[0m outdir')
	parser.add_argument('-n', metavar='NUMSAMPLES', type=int, default=1, help='How many samples to generate')
	parser.add_argument('-o', '--outdir', help='Where to place the sequences')
	parser.add_argument('--seed', help='Seed to use')
	parser.add_argument('--replacement', action='store_true', help='Draw samples with replacement')

	args = parser.parse_args()

	if args.seed: seed = int(args.seed, 16)
	else: seed = random.randint(0, 0xffffffffffffffff)

	#random.seed(seed)

	main(args.infile, samplesize=args.c, samples=args.n, outdir=args.outdir, initseed=seed, replacement=args.replacement, merge=args.merge)
