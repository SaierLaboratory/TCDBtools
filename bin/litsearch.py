#!/usr/bin/env python
from __future__ import print_function, unicode_literals

import subprocess
import argparse
import shlex
import os
import tempfile
import re
from xml.etree.ElementTree import ElementTree as ET
import urllib
import time
import sys
import io

DELAY = 0.4

def blast(blastargs, blastkwargs):
	cmd = ['blastp']
	cmd.extend([str(x) for x in blastargs])
	for k in blastkwargs: cmd.extend([str(k), str(blastkwargs[k])])
	out = subprocess.check_output(cmd)
	return out

def info(*things): print('[INFO]', *things, file=sys.stderr)

def extract_accession(teststr):
	#1. Refseq, with its distinct underscore
	refseq = re.findall('[ANYXW]P_[0-9]+(?:\.[0-9]+)', teststr)
	if refseq: 
		#print('refseq', refseq[0])
		return refseq[0]

	#2. Uniprot, a fairly common accession on TCBD
	uniprot = re.findall('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}', teststr)
	if uniprot: 
		#print('uniprot', uniprot[0])
		return uniprot[0]

	#3. Genbank
	genbank = re.findall('[A-Z]{3}[0-9]{5,7}(?:\.[0-9]+)', teststr)
	if genbank:
		return genbank[0]

	#4. PDB, liable to *bleed if detected before Uniprot accessions owing to its shorter length
	pdb = re.findall('[0-9][A-Z][0-9A-Z]{2}', teststr)
	if pdb:
		return pdb[0]

	#5. GIs, retired in theory since 2017
	gi = re.findall('[0-9]{7,}', teststr)
	if gi:
		return gi[0]

	#99. Unknown, just send it back? Raise an exception?
	#print('unidentifiable', teststr)
	return teststr

def get_subjects(tblstr):
	acclist = []
	for l in tblstr.split('\n'):
		if not l.strip(): continue
		elif l.startswith('#'): continue

		sl = l.split('\t')
		sacc = sl[1]
		acclist.append(extract_accession(sacc))
	return acclist
		
def _curl(baseurl, postdata):
	f = io.BytesIO()
	p = subprocess.Popen(['curl', '-d', postdata, baseurl], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	f.write(out)
	f.seek(0)
	return f

def process_ids(idlist, batchsize=10000, stringent=False, count_siblings=False, sibling_batch = 100):
	summaries = {}

	total = 0
	n = 0
	for i in range(0, len(idlist), 1000):
		seqlistkey, webenv = epost(idlist[i:i+1000], 'protein')
		publistkeys = elink2(webenv, seqlistkey, 'protein', 'pubmed')
		reportme = set()
		for key in publistkeys: 
			retstart = 1
			nextbatch = None
			while nextbatch != {}:
				nextbatch = esummary2(webenv, key, 'pubmed', retstart=retstart, retmax=batchsize)
				if not nextbatch: break
				reportme.update(nextbatch)
				summaries.update(nextbatch)
				retstart += batchsize
			if stringent: break

		n += 1
		info('Batch {} of {}: Added {} summaries'.format(n, (len(idlist)//1000)+1, len(reportme)))
		total += len(reportme)
	info('Finished retrieving {} summaries'.format(total))

	if count_siblings:
		info('Now counting linked proteins')
		sorsum = [str(x) for x in summaries]
		for pmid in summaries:
			n = elink_count([str(pmid)], dbfrom='pubmed', dbto='protein')
			summaries[pmid]['siblings'] = n

	return summaries

def elink2(webenv, querykey, dbfrom, dbto):
	baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
	postdata = {
		'WebEnv':webenv,
		'dbfrom':dbfrom, 
		'db':dbto,
		'query_key':querykey,
		'cmd':'neighbor_history'
	}
	#f = urllib.urlopen(baseurl, urllib.urlencode(postdata))
	f = _curl(baseurl, urllib.urlencode(postdata))

	et = ET()
	elinkresult = et.parse(f)
	query_keys = []
	for linkset in elinkresult: 
		for thing in linkset: 
			for subthing in thing:
				if subthing.tag.endswith('QueryKey'): query_keys.append(subthing.text)
	time.sleep(DELAY)
	return query_keys
	
def esummary2(webenv, query_key, db, retstart, retmax):
	baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
	postdata = {'db': db, 'WebEnv': webenv, 'query_key': query_key}
	postdata['retstart'] = retstart
	postdata['retmax'] = retmax

	#f = urllib.urlopen(baseurl, urllib.urlencode(postdata))
	f = _curl(baseurl, urllib.urlencode(postdata))

	et = ET()
	esumresult = et.parse(f)
	summaries = {}
	pmid = None
	current = None
	strkeys = set(('PubDate', 'EPubDate', 'Source', 'LastAuthor', 'Title', 'Volume', 'Issue', 'Pages', 'NlmUniqueID', 'ISSN', 'ESSN', 'RecordStatus', 'PubStatus', 'DOI', 'HasAbstract', 'PmcRefCount', 'FullJournalName', 'ELocationID', 'SO'))
	for docsum in esumresult: 
		for thing in docsum: 
			if thing.tag.endswith('Id'): 
				pmid = int(thing.text)
				summaries[pmid] = {}
				current = summaries[pmid]
			elif thing.tag.endswith('Item'):
				if thing.attrib['Name'] in strkeys:
					current[thing.attrib['Name']] = thing.text
				elif thing.attrib['Name'].endswith('List'): 
					current[thing.attrib['Name']] = []
					for element in thing: current[thing.attrib['Name']].append(element.text)
				elif thing.attrib['Name'] in ('ArticleIds', 'History'):
					current[thing.attrib['Name']] = {}
					for element in thing: current[thing.attrib['Name']][element.attrib['Name']] = element.text
	time.sleep(DELAY)
	return summaries

def efetch2(webenv, query_key, db, **kwargs):
	baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
	postdata = {'db':db, 'WebEnv':webenv, 'query_key':query_key}
	postdata.update(kwargs)
	f = _curl(baseurl, urllib.urlencode(postdata))

	et = ET()
	#efetchresult = et.parse(f)
	efetchresult = f.read()
	time.sleep(DELAY)
	return efetchresult

	
def epost(idlist, db):
	baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
	postdata = {'db': db}
	postdata['id'] = ','.join(list(idlist))

	#f = urllib.urlopen(baseurl, urllib.urlencode(postdata))
	f = _curl(baseurl, urllib.urlencode(postdata))

	et = ET()
	epostresult = et.parse(f)
	querykey, webenv = None, None
	for thing in epostresult: 
		if thing.tag.endswith('QueryKey'): querykey = thing.text
		elif thing.tag.endswith('WebEnv'): webenv = thing.text
	time.sleep(DELAY)
	return querykey, webenv

def elink_count(idlist, dbfrom, dbto):
	baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
	postdata = {'dbfrom': dbfrom, 'db': dbto}
	postdata['id'] = ','.join(list(idlist))

	f = _curl(baseurl, urllib.urlencode(postdata))

	et = ET()
	elinkresult = et.parse(f)

	n = 0
	for elinkset in elinkresult:
		for thing in elinkset:
			if thing.tag == 'LinkSetDb':
				for subthing in thing:
					if subthing.tag == 'Link': n += 1

	time.sleep(DELAY)
	return n

def elink(idlist, dbfrom, dbto):
	baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
	postdata = {'dbfrom': dbfrom, 'db': dbto}
	postdata['id'] = ','.join(list(idlist))

	#f = urllib.urlopen(baseurl, urllib.urlencode(postdata))
	f = _curl(baseurl, urllib.urlencode(postdata))

	et = ET()
	elinkresult = et.parse(f)
	linklist = set()
	for linkset in elinkresult: 
		for thing in linkset: 
			for subthing in thing: 
				if subthing.tag == 'Link': 
					for linkid in subthing: 
						linklist.add(linkid.text)
	time.sleep(DELAY)
	return linklist

	#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=pubmed&id=15718680,157427902,119703751
	#this delay will ensure this script doesn't run afoul of NCBI's rate limits
	time.sleep(0.35)

def esummary(idlist, db, start=0, batchsize=100):
	#specifically used for Pubmed and Pubmed Central
	#will probably break with other databases
	if type(idlist) not in (list, tuple): sorted_ids = sorted(idlist)
	else: sorted_ids = idlist

	baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
	postdata = {'db':db}

	summaries = {}

	for i in range(start, start+batchsize, batchsize):
		postdata.update({'id':','.join(sorted_ids[i:i+batchsize])})
		#f = urllib.urlopen(baseurl, data=urllib.urlencode(postdata))
		f = _curl(baseurl, urllib.urlencode(postdata))
		et = ET()
		esumresult = et.parse(f)
		curid = None
		curobj = {}
		for docsum in esumresult:
			for prop in docsum:
				if prop.tag == 'Id': 
					if curid: summaries[curid] = curobj
					curid = prop.text
					curobj = {}
				elif prop.tag == 'Item': 
					if prop.attrib['Name'] == 'AuthorList':
						authorlist = [author.text for author in prop]
						curobj['AuthorList'] = authorlist
					else: curobj[prop.attrib['Name']] = prop.text
		if curid: summaries[curid] = curobj

		time.sleep(0.35)
	return summaries


def fmtline(dbid, obj):
	out = ''
	out += '{}'.format(dbid)
	#out += '\t{}'.format(obj['PubDate'])
	out += '\t{}'.format(obj['PubDate'][:4])
	out += '\t{}'.format(obj['PmcRefCount'])
	if 'siblings' in obj: out += '\t{}'.format(obj['siblings'])
	out += '\t{}'.format(obj['Title'])
	out += '\t' + ', '.join(obj['AuthorList'])
	#if obj['AuthorList'].strip():
	#	out += '{}'.format(obj['AuthorList'])
	#	out += ', '
	#out += obj['LastAuthor']
	return out

def main(blastargs, blastkwargs, email, outfile='/dev/stdout', batchsize=5000, loadblast=None, stringent=False, count_siblings=False):

	if loadblast:
		info('Loading BLAST table "{}"'.format(blastkwargs['-query']))
		out = open(blastkwargs['-query']).read()
		newout = None
		for l in out.split('\n'):
			if l.startswith('#'): continue
			elif not l.strip(): continue
			elif '\t' not in l[:20]:
				if newout is None: newout = ''
				newout += '\t{}\n'.format(l)
		out = out if newout is None else newout
				
	else:
		info('Blasting...')
		out = blast(blastargs, blastkwargs)

	allsubjects = set()

	for sacc in get_subjects(out): allsubjects.add(sacc)
	info('Found {} unique hits'.format(len(allsubjects)))
	allsubjects = sorted(allsubjects)
	megabatch = 5000

	summaries = process_ids(allsubjects, stringent=stringent, count_siblings=count_siblings)
	info('Linked to {} articles'.format(len(summaries)))
	with open(outfile, 'w') as f:
		f.write('#Accession\tYear\tCited by')
		if count_siblings: f.write('\tLinks to')
		f.write('\tTitle\tAuthors\n')
		for pmid in sorted(summaries):
			out = fmtline(pmid, summaries[pmid])
			f.write(out.encode('utf-8'))
			f.write('\n')

	info('Done!')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--db',  default='nr', help='Database to use (default: nr)')

	if 'ENTREZ_EMAIL' in os.environ: parser.add_argument('--email', default=os.environ['ENTREZ_EMAIL'], help='Email address for Entrez queries (default: $ENTREZ_EMAIL = {})'.format(os.environ['ENTREZ_EMAIL']))
	else: parser.add_argument('--email', help='Email address for Entrez queries. Set $ENTREZ_EMAIL to avoid having to specify this each run.')

	parser.add_argument('--evalue', default=1e-10, type=float, help='E-value cutoff (default: 1e-10)')
	parser.add_argument('--local', action='store_true', help='BLAST locally. (Pubmed searches will be remote regardless)')
	parser.add_argument('--blastargs', help='Arguments to pass to BLASTP. Make sure this is enclosed in quotes.')
	parser.add_argument('--loadblast', action='store_true', help='Load a BLAST results table instead of BLASTing (or anything tab-separated with #-commented lines and subject accessions in the second column)')
	parser.add_argument('--stringent', action='store_true', help='Use only the best articles during the Entrez Elink step. Greatly improves performance and precision at the expense of sensitivity')
	parser.add_argument('--count-siblings', action='store_true', help='Count proteins linked to each PMID (!!VERY SLOW!!)')
	parser.add_argument('infile', help='File containing sequences of interest. Use "asis:" to enter raw sequences and "stdin" or "-" to read from standard input')
	parser.add_argument('outfile', nargs='?', default='stdout', help='Where to write the collected Pubmed and Pubmed Central results')
	
	args = parser.parse_args()

	if not args.email: raise TypeError('Please specify an email address for Entrez queries (--email or set $ENTREZ_EMAIL)')

	if (args.infile) in ('-', 'stdin'): 
		info('Now reading from stdin')
		infile = '/dev/stdin'
	elif args.infile.startswith('asis:'):
		f = tempfile.NamedTemporaryFile()
		f.write(args.infile[5:])
		f.flush()	
	else: infile = args.infile

	if args.outfile in ('-', 'stdout'):
		outfile = '/dev/stdout'
	else: outfile = args.outfile


	blastargs = []
	blastkwargs = {}


	blastkwargs.update({'-query':infile})
	blastkwargs.update({'-evalue':args.evalue})
	blastkwargs.update({'-db':args.db})
	blastkwargs.update({'-outfmt':6})
	if not args.local: blastargs.append('-remote')

	#blastargs.extend(shlex.split(args.blastargs))


	if args.blastargs: reqargs = shlex.split(args.blastargs)
	else: reqargs = []
	i = 0
	while reqargs:
		k = reqargs.pop(0)
		if reqargs[0].startswith('-'): blastargs.append(k)
		else: blastkwargs.update({k:reqargs.pop(0)})

	main(blastargs, blastkwargs, args.email, outfile, loadblast=args.loadblast, stringent=args.stringent, count_siblings=args.count_siblings)

