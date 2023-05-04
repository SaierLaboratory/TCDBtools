#!/usr/bin/env python2

from __future__ import print_function

import sys, subprocess, os, time

def usage():
	print("Usage: {0} tcid genomeid".format(sys.argv[0]))
	exit(1)

def makedb(fn, outdir):
	title = os.path.splitext(os.path.basename(fn))[0]
	outpath = outdir + '/' + title
	title += ' {}'.format(time.strftime('%Y-%m-%d %H:%M:%S'))
	try: subprocess.call(['makeblastdb', '-dbtype', 'prot', '-parse_seqids', '-hash_index', '-title', title, '-in', fn, '-out', outpath])
	except subprocess.CalledProcessError:
		print('Could not run makeblastdb')
		return 1
	#makeblastdb -dbtype prot -parse_seqids -hash_index -title test -in /ResearchData/Users/nuo/GenomeAnalysis/GCA_001029695.1_ASM102969v1_protein.faa -out nuo
def extractfamily(outdir):
	try: out = subprocess.check_output(['extractFamily.pl', '-i', 'all', '-o', outdir, '-f', 'blast'])
	except subprocess.CalledProcessError: 
		print('Could not run extractFamily.pl')
		return 1

def run_pfam(tcid, afmid, outdir, db):
	fn = '{0}/family-{1}.faa'.format(outdir, tcid)

	#if os.path.isfile(fn) and os.path.getsize(fn): pass
	#else:
	#	try: out = subprocess.check_output(['extractFamily.pl', '-i', tcid, '-o', outdir + '/blast', '-f', 'blast'])
	#	except subprocess.CalledProcessError: 
	#		print('Could not run extractFamily.pl')
	#		return 1

	#try: out = subprocess.check_output(['getSequence', afmid])
	out = ''
	try: out += subprocess.check_output(['blastdbcmd', '-db', db, '-entry', afmid, '-dbtype', 'prot'])
	except subprocess.CalledProcessError: 
		print('Could not run getSequence for {} in {}'.format(afmid, db))
		return 1
	out += '\n'
	try: out += subprocess.check_output(['blastdbcmd', '-db', 'tcdb', '-entry', tcid, '-dbtype', 'prot'])
	except subprocess.CalledProcessError: 
		print('Could not run getSequence for {} in tcdb'.format(tcid))
		return 1
	with open(fn, 'a') as f: f.write(out)

	#try: out = subprocess.check_output(['run_pfam', fn, '{0}/family-{1}.pfam'.format(outdir, tcid)])
	#except subprocess.CalledProcessError: 
	#	print('Could not run run_pfam')
	#	return 1

	try: 
		pfamout = subprocess.check_output(['hmmscan', '--cpu', '4', '--noali', '--cut_ga', '-o', '/dev/null', '--domtblout', '{}/family-{}.pfam'.format(outdir, tcid), '/ResearchData/pfam/pfamdb/Pfam-A.hmm', fn])
		
		#with open(outdir + '/{}.pfam'.format(tcid), 'w') as f: f.write(pfamout)
		print('Your results are in: {0}/family-{1}.pfam'.format(outdir, tcid))
	except subprocess.CalledProcessError:
		print('Could not run hmmscan for {} and {}'.format(tcid, afmid))


def interpret(line, outdir='.', db='tcdb'):
	if not line.strip(): return 0
	elif line.lstrip().startswith('#'): return 0
	elif len(line.split('\t')) >= 3: 
		#default: genome_acc ??? TCDB_acc
		tcid = line.strip().split()[2] + '-' + line.strip().split()[1]
		afmid = line.strip().split()[0]
		run_pfam(tcid, afmid, outdir=outdir, db=db)
		return 0
	else: 
		print("Ask Kevin to implement 2-arg handling (it'll take only a minute or two)")
		return 1

def main(infile, genomefn, outdir='.'):
	makedb(genomefn, outdir + '/blast')
	extractfamily(outdir + '/blast')

	for fn in os.listdir(outdir):
		if fn.startswith('family-') and fn.endswith('.faa'): os.remove(outdir + '/' + fn)

	db = os.path.splitext(os.path.basename(genomefn))[0]

	#while 1:
		#interpret(raw_input('> '), outdir=outdir, db=db)
	with open(infile) as f:
		for l in f: interpret(l, outdir=outdir, db=db)

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()

	#parser.add_argument('-c', action='store_true', help='continuous mode')
	parser.add_argument('infile', help='genome hits table (generally something.tsv) ')
	parser.add_argument('-d', default='tcdb', help='proteome fasta (generally something.faa)')
	parser.add_argument('-o', default='afhbot_out', help='directory to dump files in (default: afhbot_out)')

	args = parser.parse_args()

	if 'BLASTDB' in os.environ:
		os.environ['BLASTDB'] = args.o + '/blast:' + os.environ['BLASTDB']
	else:
		os.environ['BLASTDB'] = args.o + '/blast'

	if not os.path.isdir(args.o): os.mkdir(args.o)
	if not os.path.isdir(args.o + '/blast'): os.mkdir(args.o + '/blast')

	main(infile=args.infile, genomefn=args.d, outdir=args.o)

##turbo mode
#[ $1 ] || usage
#[ $2 ] || usage
#[ $3 ] || usage
#
#extractFamily.pl -i $3 -o .
#getSequence $1 >> family-$3.faa
#run_pfam family-$3.faa family-$3.pfam
#
#echo "Your results are in: family-$3.pfam"
