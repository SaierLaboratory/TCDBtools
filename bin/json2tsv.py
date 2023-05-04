#!/usr/bin/env python
from __future__ import print_function
import json, argparse, sys
import os

def pretty(obj, k):
	if 'cov' in k: return '{:0.2%}'.format(obj[k])
	elif 'present' in k: return '{}-{}'.format(obj[k][0][0], obj[k][0][1])
	else: return str(obj[k])

def dump_file(fin, fout, initial=True, blacklist=None, whitelist=None):
	header = ''
	blacklist = set(['distances', 'qaligned', 'saligned', 'queryfn', 'subjectfn']) if blacklist is None else blacklist
	#length	matrix	qhel	qpresent	quality	query	rmsd	shel	spresensubject
	whitelist = ['name', 'query', 'qchain', 'qhel', 'subject', 'schain', 'shel', 'rmsd', 'quality', 'length', 'tmcov', 'fullcov', 'qtmlen', 'qlen', 'stmlen', 'slen', 'qpresent', 'spresent'] if whitelist is None else whitelist
	for l in fin:
		if not l.strip(): continue
		elif l.startswith('#'): continue
		name, rawjson = l.split('\t')
		query, qchain, qhel, vs, subject, schain, shel = name.split('_')
		nameobj = {}
		nameobj['name'] = name
		nameobj['query'] = query
		nameobj['qchain'] = qchain
		nameobj['qhel'] = qhel
		nameobj['subject'] = subject
		nameobj['schain'] = schain
		nameobj['shel'] = shel
		obj = json.loads(rawjson)

		if not header:
			#header += '#name\t'
			#header = '#name\tquery\tqchain\tqtms\tsubject\tschain\tstms\t'
			header = '#'
			for k in whitelist:
				if k in blacklist: continue
				header += k + '\t'
			fout.write(header.strip() + '\n')
		s = ''
		#s += '{}\t'.format(name)
		#s += '{}\t{}\t{}\t'.format(query, qchain, qhel)
		#s += '{}\t{}\t{}\t'.format(subject, schain, shel)
		for k in whitelist:
			if k in blacklist: continue
			elif k in nameobj: s += nameobj[k] + '\t'
			else:
				try: s += pretty(obj, k) + '\t'
				except KeyError: s += '\t'
		fout.write(s.strip() + '\n')
		

def main(infile, outfile, blacklist=None, whitelist=None):
	if isinstance(outfile, str): fout = open(outfile, 'w')
	else: fout = outfile

	if os.path.isfile(infile):
		fin = open(infile)
		dump_file(fin, fout, blacklist=blacklist, whitelist=whitelist)
	elif os.path.isdir(infile):
		initial = True
		for d2subdir in os.listdir(infile):
			fn = '{}/{}/tmalignments/sp_all.tsv'.format(infile, d2subdir)
			if os.path.isfile(fn):
				fin = open(fn)
				dump_file(fin, fout, initial=initial, blacklist=blacklist, whitelist=whitelist)
				initial = False

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', nargs='?', default='/dev/stdin', help='input file or Deuterocol2 directory (default:stdin)')
	parser.add_argument('outfile', nargs='?', default='/dev/stdout', help='output file (default:stdout)')
	parser.add_argument('--blacklist', nargs='*', help='Which fields to exclude')
	parser.add_argument('--whitelist', nargs='*', help='Which fields to include')

	args = parser.parse_args()
	if args.infile == '/dev/stdin': 
		print('Now reading input from stdin...', file=sys.stderr)

	main(args.infile, args.outfile, blacklist=args.blacklist, whitelist=args.whitelist)
