#!/usr/bin/env python

from Bio import SeqIO
import argparse
import re

def main(pattern, fh, seqmode=False, invert=False, single=False):
	recording = False

	results = []
	for record in SeqIO.parse(fh, 'fasta'):
		addme = False
		if seqmode:
			if invert ^ bool(re.search(pattern, str(record.seq))): addme = True
		else:
			if record.name and not record.description: searchme = record.name
			elif record.description and not record.name: searchme = record.description
			elif record.name and record.description: searchme = record.name
			else: searchme = ''

			if invert ^ (bool(searchme) and bool(re.search(pattern, searchme))): addme = True
			#elif invert ^ (bool(record.description) and bool(re.search(pattern, record.description))): addme = True
		if addme: 
			results.append(record)
			if single: break
	return results

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-c', action='store_true', help='Count matches')
	parser.add_argument('-f', action='store_true', help='Treat PATTERN as a filename containing newline-separated patterns')
	parser.add_argument('-o', action='store_true', help='Retrieve only the accessions (and not the sequences)')
	parser.add_argument('-s', action='store_true', help='Sequence grep mode, i.e. grep in sequences rather than in headers')
	parser.add_argument('-v', action='store_true', help='Inverse match')
	parser.add_argument('--single', action='store_true', help='Stop at first match')
	parser.add_argument('pattern', help='Pattern to grep for')
	parser.add_argument('infile', nargs='*', default=['/dev/stdin'], help='Files to read in (default: stdin)')

	args = parser.parse_args()

	if not args.f: pattern = args.pattern
	else: 
		patternlist = []
		with open(args.pattern) as fh:
			for l in fh: patternlist.append(l.replace('\n', ''))
		pattern = '|'.join(patternlist)

	for fn in args.infile:

		with open(fn) as f: results = main(pattern, f, seqmode=args.s, invert=args.v, single=args.single)
		out = ''

		if args.c: print(len(results))
		else:
			for seq in results: 
				if not args.o:
					if seq.description: print('>{}\n{}'.format(seq.description, seq.seq))
					elif seq.name: print('>{}\n{}'.format(seq.name, seq.seq))
					else: print('>\n{}'.format(seq.seq))
				else:
					if seq.description: print('{}'.format(seq.description))
					elif seq.name: print('{}'.format(seq.name))
					else: print('')
