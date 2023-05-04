#!/usr/bin/env python

import argparse
import os
import sys
import entrdbcmd
from Bio import SeqIO
import subprocess
import shlex
import re
import intrafam
import numpy as np

def info(*things):
	print('[INFO]', *things, file=sys.stderr)
def warn(*things):
	print('[WARNING]', *things, file=sys.stderr)
def error(*things):
	print('[ERROR]', *things, file=sys.stderr)
	exit(1)

def collect_blast_results(fh, expect_cutoff, column=10):
	passing_pairs = set()
	for l in fh:
		if not l.strip(): continue
		elif l.startswith('#'): continue

		sl = l.replace('\n', '').split('\t')

		query, subject = sl[:2]

		try: expect = float(sl[column])
		except IndexError: error('Not enough columns in BLAST table (requested {}, found {})'.format(column, len(sl)))
		except ValueError: error('Invalid e-value: {}'.format(sl[column]))

		if expect > expect_cutoff: continue
		else: passing_pairs.add((query, subject))

	return passing_pairs

def extract_subjects(pairlist):
	subjects = set()
	for pair in pairlist: subjects.add(pair[1])
	return subjects

def get_yn(prompt, default='y'):
	s = prompt
	if default == 'y': s += '[Y/n] '
	else: s += '[y/N] '

	resp = input(prompt)

	if resp.strip().lower().startswith('y'): return True
	elif resp.strip().lower().startswith('n'): return False
	elif default == 'y': return True
	else: return False

def get_seqlengths(fh):
	seqlengths = {}
	for seqobj in SeqIO.parse(fh, 'fasta'):
		seqlengths[seqobj.name] = len(seqobj.seq)
	fh.seek(0)
	return seqlengths

def mass_align(fn1, fn2, align_tool='ssearch36', args=None, outdir=None, expect_cutoff=1e-3):
	if align_tool.startswith('ssearch'):
		cmd = [align_tool, '-m', '8', '-E', str(expect_cutoff)]
		if args: cmd.extend(args)
		cmd.extend([fn1, fn2])
		out = subprocess.check_output(cmd).decode('utf-8')
	elif align_tool.startswith('blastp'):
		cmd = [align_tool, '-outfmt', '6', '-query', fn1, '-subject', fn2, '-evalue', str(expect_cutoff)]
		if args: cmd.extend(args)
		out = subprocess.check_output(cmd).decode('utf-8')

	elif align_tool.startswith('swipe'):
		subprocess.call(['makeblastdb', '-in', fn2, '-dbtype', 'prot', '-out', '{}/swipedb'.format(outdir), '-hash_index', '-parse_seqids'])
		cmd = [align_tool, '-i', fn1, '-d', '{}/swipedb'.format(outdir), '-m', '8', '-e', str(expect_cutoff)]
		if args: cmd.extend(args)
		out = subprocess.check_output(cmd).decode('utf-8')

		out = re.sub('[a-z]+\|([^\t|]+)\|', r'\1', out)
		

	else: error('Unsupported alignment tool: {}'.format(align_tool))

	return out

def cdhit(infn, outfn, ident=0.9, cov=0.8):
	cmd = ['cd-hit', '-i', infn, '-o', outfn, '-c', str(ident), '-G', '0', '-aS', str(cov)]

	subprocess.call(cmd)

	representatives = set()
	with open(outfn) as fh:
		for record in SeqIO.parse(fh, 'fasta'):
			representatives.add(record.name)
	os.rename(outfn, outfn + '.faa')
	return representatives

def parse_methodstr(methodstr):
	kwargs = {}
	args = methodstr.split(':')
	skip = 0
	for i, arg in enumerate(args):
		if skip > 0: 
			skip -= 1
			continue

		if arg in ('fixed'): 
			kwargs['scorefunction'] = intrafam.ScoreFunction.get_method(arg)
			kwargs['threshold'] = float(args[i+1])
			skip += 1

		elif arg == 'selfnorm': kwargs['normalization'] = 'selfnorm'

	if 'scorefunction' not in kwargs: kwargs['scorefunction'] = intrafam.ScoreFunction.fixed
	elif 'normalization' not in kwargs: kwargs['normalization'] = 'selfnorm'

	return kwargs

def score(obj, scorefunction):
	q, s = scorefunction(obj)

	qdelta = np.zeros(obj['qlen'])
	qdelta[obj['qspan'][0]-1:obj['qspan'][1]] += q
	sdelta = np.zeros(obj['slen'])
	sdelta[obj['sspan'][0]-1:obj['sspan'][1]] += s
	return qdelta, sdelta

def cumulate(blastfn, pairs, seqlen, expect_cutoff=1e-3, scorefunction='fixed', normalization='selfnorm', threshold=0.7, **kwargs):# method='fixed:selfnorm:0.7:grablargest:0.2'):
	if isinstance(scorefunction, str): scorefunction = intrafam.ScoreFunction.get_method(scorefunction)

	rawcounts = {}

	with open(blastfn) as fh:
		for l in fh: 
			if not l.strip(): continue
			elif l.startswith('#'): continue
			sl = l.replace('\n', '').split('\t')

			if float(sl[10]) > expect_cutoff: continue
			query = sl[0]
			subject = sl[1]
			pident = float(sl[2])
			qspan = (int(sl[6]), int(sl[7]))
			sspan = (int(sl[8]), int(sl[9]))
			expect = float(sl[10])
			bits = float(sl[11])
			qname, sname = sl[:2]
			qlen = seqlen[query]
			slen = seqlen[subject]

			obj = {'query':query, 'subject':subject,
				'pident':pident,
				'qspan':qspan, 'sspan':sspan,
				'qlen':qlen, 'slen':slen}

			if qname not in rawcounts: rawcounts[qname] = np.zeros(seqlen[qname])
			if sname not in rawcounts: rawcounts[sname] = np.zeros(seqlen[sname])

			#rawcounts[qname][qspan[0]-1:qspan[1]] += 
			qdelta, sdelta = score(obj, scorefunction)

			rawcounts[qname] += qdelta
			rawcounts[sname] += sdelta
	#TODO: extend for other normalization schemes

	if normalization == 'selfnorm':
		for seqname in rawcounts:
			if max(rawcounts[seqname]) != 0: rawcounts[seqname] /= max(rawcounts[seqname])

	return rawcounts

def grab_contigs(arr, threshold=0.7, decay=0.8, permitlength=150):
	last = None
	contigs = []
	for i, v in enumerate(arr):
		if (v >= threshold):
			if (i-1) != last: contigs.append([i, i])
			else: contigs[-1][-1] += 1
			last = i

	if not contigs: return contigs

	lengths = [span[1] - span[0] + 1 for span in contigs]

	maxlength = max(lengths)

	finalcontigs = []
	for length, span in zip(lengths, contigs):
		if length >= permitlength: 
			finalcontigs.append(span)
		elif length >= (decay * maxlength):
			finalcontigs.append(span)
	return finalcontigs

def main(infile, outdir, fxpand_expect=1e-20, prune_expect=1e-20, cdhit_ident=0.8, cdhit_cov=0.8, method='fixed:0.7:selfnorm', force=False, column=10, align_tool='ssearch36', align_tool_args=None, interactive=False, threshold=0.7, decay=0.8, permitlength=150):

	if os.path.isdir(outdir) and not force: warn('Directory `{}\' already exists'.format(outdir))
	else: os.mkdir(outdir)

	with open(infile) as fh:
		fxpand_pairs = collect_blast_results(fh, fxpand_expect, column)
		info('Found {} alignments'.format(len(fxpand_pairs)))
	subjects = extract_subjects(fxpand_pairs)
	info('Found {} subject sequences'.format(len(subjects)))

	if interactive: 
		proceed = get_yn('Continue by fetching {} sequences?'.format(len(subjects)), default='y')
		if proceed: pass
		else: exit(0)
		

	if os.path.isfile('{}/subjects.faa'.format(outdir)) and not force:
		warn('Existing subjects.faa found. Skipping efetch step. Remove it or use --force to overwrite it.')
	else:
		with open('{}/subjects.faa'.format(outdir), 'w') as fh:
			fh.write(entrdbcmd.efetch(list(subjects)))

	seqlengths = {}
	with open('{}/subjects.faa'.format(outdir)) as fh:
		seqlengths = get_seqlengths(fh)


	representatives = cdhit('{}/subjects.faa'.format(outdir), '{}/clustered'.format(outdir), ident=cdhit_ident, cov=cdhit_cov)

	
	if os.path.isfile('{}/allvall.tsv'.format(outdir)) and not force:
		warn('Existing allvall.tsv found. Skipping align step. Remove it or use --force to overwrite it.')
	else: 
		allvallstr = mass_align('{}/clustered.faa'.format(outdir), '{}/clustered.faa'.format(outdir), align_tool=align_tool, args=align_tool_args, outdir=outdir, expect_cutoff=1e-3)
		with open('{}/allvall.tsv'.format(outdir), 'w') as fh:
			fh.write(allvallstr)

	with open('{}/allvall.tsv'.format(outdir)) as fh:
		pairs = collect_blast_results(fh, prune_expect)


	contigs = {}
	overlapscores = cumulate('{}/allvall.tsv'.format(outdir), pairs=pairs, seqlen=seqlengths, expect_cutoff=prune_expect, **parse_methodstr(method))
	for seq in overlapscores:
		contigs[seq] = grab_contigs(overlapscores[seq], threshold=threshold, permitlength=permitlength)

	outfh = open('{}/contigs.faa'.format(outdir), 'w')
	with open('{}/clustered.faa'.format(outdir)) as fh:
		for record in SeqIO.parse(fh, 'fasta'):
			if record.name in contigs and len(contigs[record.name]):
				for i, span in enumerate(contigs[record.name]):
					header = '>{}_contig{}'.format(record.name, i)
					seq = record.seq[span[0]:span[1]+1]
					outfh.write('{}\n{}\n'.format(header, seq))

	info('Successfully wrote contigs to {}/contigs.faa'.format(outdir))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	#BLAST filtering params
	parser.add_argument('-k', type=int, default=11, help='1-indexed to match `cut\'/`sort\'. Which column contains the e-value (if not the default, 11. Use -k5 for famXpander-written tables)')

	parser.add_argument('-e1', '--fxpand-expect', type=float, default=1e-20, help='Expect threshold for alignments in the input (default: 1e-20)')

	parser.add_argument('-e2', '--prune-expect', type=float, default=1e-20, help='Expect threshold for alignments from the all-vs-all step (default: 1e-20)')

	parser.add_argument('--align-tool', default='ssearch36', help='Program to use for the all-vs-all alignments (default: ssearch36)')
	parser.add_argument('--align-tool-args', default='', help='Arguments to pass to align tool')

	parser.add_argument('--threshold', type=float, default=0.7, help='Relative occupancy threshold for contiguous segments')
	parser.add_argument('--decay', type=float, default=0.8, help='Minimum relative length of accepted contigs')
	parser.add_argument('--permit-length', type=int, default=150, help='Permit all contigs as long as or longer than this (default: 150)')

	#CD-hit params
	parser.add_argument('-i', type=float, default=0.8, help='Identity cutoff for CD-hit')

	parser.add_argument('-c', type=float, default=0.8, help='Coverage (on smaller seq.) for CD-hit')

	#block selection params
	parser.add_argument('-m', '--method', default='fixed:0.7', help='Block selection method string (default: fixed:0.7:selfnorm)')

	#misc. args
	parser.add_argument('--force', action='store_true', help='Force this program to overwrite existing data')
	parser.add_argument('--interactive', action='store_true', help='Prompt at every significant step')

	#positionals
	parser.add_argument('infile', help='BLAST-like TSV')
	parser.add_argument('-o', required=True, help='Output directory')

	args = parser.parse_args()

	main(args.infile, args.o, 
		fxpand_expect=args.fxpand_expect, prune_expect=args.prune_expect, column=args.k-1,
		align_tool=args.align_tool, align_tool_args=shlex.split(args.align_tool_args),
		cdhit_ident=args.i, cdhit_cov=args.c,
		method=args.method, force=args.force, 
		threshold=args.threshold, decay=args.decay, permitlength=args.permit_length,
		interactive=False)
