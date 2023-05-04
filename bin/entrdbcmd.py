#!/usr/bin/env python

import argparse
import os
import subprocess

EMAIL = ''

def guess_database(dbname, stopcrashingplease=False):
	if dbname == 'nr': return 'protein'
	elif dbname == 'nt': return 'nucleotide'
	elif stopcrashingplease: return 'protein'
	else: raise ValueError('Invalid db name: `{}\''.format(dbname))

def _efetch(acclist, db='protein'):
	postdata = ''
	postdata += 'db={}'.format(db)
	postdata += '&retmode=text'
	postdata += '&rettype=fasta'
	postdata += '&id='
	postdata += ','.join(acclist)
	postdata = postdata.encode('utf-8')

	url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

	p = subprocess.Popen(['curl', url, '-d', '@-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)

	out, err = p.communicate(input=postdata)
	out = out.decode('utf-8')
	out += '\n'
	return out


def efetch(acclist, db='protein', batchsize=1000):
	sequences = ''
	for i in range(0, len(acclist), batchsize):
		sequences += _efetch(acclist[i:i+batchsize], db)

	return sequences

def load_acclist(fh):
	acclist = []
	for l in fh: acclist.append(l.strip())
	return acclist

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-db', default='nr')

	if 'ENTREZ_EMAIL' in os.environ: parser.add_argument('-email', default=os.environ['ENTREZ_EMAIL'], help='Entrez email. Will read from $ENTREZ_EMAIL if set. (optional, default: $ENTREZ_EMAIL = {})'.format(os.environ['ENTREZ_EMAIL']))
	else: parser.add_argument('-email', required=True, help='Entrez email. Will read from $ENTREZ_EMAIL if set. (required)')


	parser.add_argument('-target_only', action='store_true', help='Does nothing. Included to allow drop-in replacement')

	entryargs = parser.add_mutually_exclusive_group(required=True)
	entryargs.add_argument('-entry', help='\033[1mCOMMA\033[0m-separated list of accessions to retrieve')
	entryargs.add_argument('-entry_batch', help='Load a \033[1mNEWLINE\033[0m-separated list of accessions')

	parser.add_argument('-out', default='/dev/stdout', help='Where to write the sequences')

	args = parser.parse_args()


	if args.entry is None: acclist = load_acclist(open(args.entry_batch))
	else: acclist = args.entry.split(',')

	db = guess_database(args.db)
	EMAIL = args.email
	sequences = efetch(acclist, db=db)

	print(sequences)
	#with open(args.out, 'w') as fh
