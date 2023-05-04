#!/usr/bin/env python

'''
 * File name: searchMissingComponent.py
 * Author: Yichi Zhang
 * Date: 2018/7/22
 * Description: This program is provided to Dr. Saier's Lab written in Python2.7 and biopython package.
 *		It suppose to let user input keywords of a protein and retreive relavant sequences from NCBI to make a blast database. It will perform blastp to form the blast result. The table blast result will be filtered and form the final html result with corresponding alignment hydropathy graphs. Pls upgrade the local biopython and pandas. Install certifi package
 * Parameters: keyword - the target protein description
 *	       genome - complete protein sequences of a genome
 * 	       evalue - evalue threashold for blast
 *	       coverage - the percent coverage of the smaller protein
 *	       out - name of the output file
 *	       o - overwrite option
 *             l - enable the sequences to be retrieved locally
'''

import os
from StringIO import StringIO
import sys
import getpass
import subprocess
import argparse
import ssl
from Bio import Entrez, SeqIO
import pandas as pd
import requests
import certifi
from hmmscanParser import hmmParser
import pprint
import bz2
import gzip
import tempfile
from Bio.Blast import NCBIXML
from formatQuodDomain import Domain_string

__author__ = 'Yichi Zhang'
__date__ = 'August 4, 2018'

# create a list containing 7 lists, each of 5 elements, all set to 0
# colors = [[0 for x in range(5)] for y in range(7)]
colors = ['red','teal','blue','yellow','purple','maroon','black','pink','green','cyan','beige','lavendar','orange','grey',
	  'olive','navy','magenta','coral','lime','plum','chocolate','indigo','brown']
'''
 * Function name: retrieve_seq( keywords, out_directory )
 * Description: search and get desired sequences from NCBI and form a txt file. It should be used in a loop to process mutilple sequences
 * Parameter: keyword must be a string. The name of the output directory
 * Return type: true or false to indicate if the search is successful or not
'''
def retrieve_seq( keywords, out_directory, flag, overwrite, redundancy ):
	cwd = os.getcwd()
	if not overwrite and os.path.isfile( '{c}/{o}/NCBI_seq.fasta'.format( c=cwd, o=out_directory ) ):
		return True
	# bypass ssl different than the lab computer
	if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
		ssl._create_default_https_context = ssl._create_unverified_context
	Entrez.email = "yiz370@ucsd.edu"
	handle = Entrez.esearch( db = "protein", term = keywords, retmax = '1000000' )
	record = Entrez.read( handle )
	# check if the results are greater than 0
	if record["Count"] == 0:
		return False
	else:
		# get the current working directory
		outpath = '{c}/{o}'.format( c=cwd, o=out_directory )
		if not os.path.isdir(outpath):
			os.system( 'mkdir {p}'.format( p=outpath ) )
		# get the first 1000,000 sequence id
		idlist = record["IdList"]
		if flag == True:
			# get the accession list
			list_to_file('{o}/NCBI_acc_list.txt'.format(o=outpath), idlist)
			os.system( 'blastdbcmd -db nr -entry_batch {o}/NCBI_acc_list.txt -target_only -out {o}/NCBI_seq.fasta'.format( o=outpath ) )
		else:
			handle2 = Entrez.efetch( db = "protein", id = idlist, rettype = "fasta", retmode = "text" )
			with open( outpath+'/NCBI_seq_raw.fasta', 'w' ) as outfile:
				for seq_record in SeqIO.parse(handle2,"fasta"):
					SeqIO.write(seq_record, outfile, "fasta")
			outfile.closed
			# check if the number of sequence is complete or not. If the percentage < 90, download locally and compare with the exsisting file, keep the one with more sequences
			with open( outpath+'/NCBI_seq_raw.fasta', 'r' ) as infile:
				# count the number of '>' in the file
				text = infile.read().strip()
				freq = text.count('>')
			infile.closed
			# compare freq os NCBI and the freq of the local file
			if float(freq)/float(record['Count']) < 0.9:
				# download locally and compare with the Internet version
				print "Not enough sequences retrieved from the website! Try to retrieve locally."
				list_to_file('{o}/NCBI_acc_list.txt'.format(o=outpath), idlist)
				os.system('blastdbcmd -db nr -entry_batch {o}/NCBI_acc_list.txt -target_only -out {o}/NCBI_seq_raw_local.fasta'.format( o=outpath ))
				with open( outpath+'/NCBI_seq_raw_local.fasta'.format( o=outpath ) ) as local_file:
					text2 = local_file.read().strip()
					freq2 = text2.count('>')
				local_file.closed
				if freq > freq2:
					os.system('rm {o}/NCBI_seq_raw_local.fasta'.format( o=outpath ))
				else:
					os.system('rm {o}/NCBI_seq_raw.fasta'.format( o=outpath ))
					os.system('mv {o}/NCBI_seq_raw_local.fasta {o}/NCBI_seq_raw.fasta'.format( o=outpath ))
		# remove redundancy
		os.system( 'cd-hit -i {o}/NCBI_seq_raw.fasta -o {o}/NCBI_seq.fasta -c {r}'.format( o=outpath, r=redundancy ) ) 
		return True

'''
 * Function name: create_db( genome, out_directory )
 * Description: compare the number of retrieved sequences with the number of local genome proteins and create a blast database for the one with more sequences
 * Parameter: genome to provide the genome file; out_directory indicates the ouput folder
 * Return type: True - GENOME DB; FALSE - NCBI DB
'''
def create_db( genome, out_directory ):
	# count the number of proteins in the genome file
	cwd = os.getcwd()
	path2 = '{c}/{o}/NCBI_seq.fasta'.format( c=cwd, o=out_directory )
	process = subprocess.Popen(genome+'grep -c ">"', stdout=subprocess.PIPE, shell=True)
	local_count = process.stdout.read().split()
	# count the number of proteins in the NCBI file
	NCBI = open( path2, 'r' )
	NCBI_list = NCBI.read()
	NCBI_count = NCBI_list.count('>')
	NCBI.close()
	#compare and create the local database
	#os.system( 'mkdir {c}/{o}/blastdb' )
	if local_count > NCBI_count:
		cmdstring = '{g}makeblastdb -in - -dbtype "prot" -title "testGenome" -out {c}/{o}/blastdb/genome_blastdb -hash_index -parse_seqids;makeblastdb -in {i2} -dbtype "prot" -out {c}/{o}/blastdb/NCBI_blastdb -hash_index -parse_seqids'.format( g=genome, c=cwd, o=out_directory, i2=path2 )
		result = 'genome'
	else:
		cmdstring = 'makeblastdb -in {i} -dbtype "prot" -out {c}/{o}/blastdb/NCBI_blastdb -hash_index -parse_seqids;{g}makeblastdb -in - -dbtype "prot" -title "testGenome" -out {c}/{o}/blastdb/genome_blastdb -hash_index -parse_seqids'.format( i=path2, c=cwd, o=out_directory, g=genome )
		result = 'NCBI'
	os.system( cmdstring )
	return result

'''
 * Function name: get_output( query, dbname, evalue, coverage, e_filter, domain_coverage, out, overwrite)
 * Description: this function will get the original blast output and parse to pandas to filter once. this result will ran pfam and filter again. Later it will form the hydropathy plot and get the final html result
 * Parameter: flags used for blast and filter
 * Return type: None
'''
def get_output( display ,dcov, overwrite, query, dbname, evalue, coverage, which, e_filter, out, alignment_length ):
	# form the blast result first
	if 'genome' in dbname:
		cmdstring = 'blastp -query {q} -db {d} -evalue {e} -max_hsps 1 -out {o}/result.xml -outfmt "5"'.format( q=query, d=dbname, e=evalue, o=out)
		col_list = ['sacc','qacc', 'slen', 'qlen', 'length','gaps', 'sstart', 'send', 'qstart', 'qend', 'evalue', 'score', 'pident', 'sseq', 'match','qseq']
	else:
		cmdstring = '{q}blastp -query - -db {d} -evalue {e} -max_hsps 1 -out {o}/result.xml -outfmt "5"'.format( q=query, d=dbname, e=evalue, o=out)
		col_list = ['qacc','sacc', 'qlen', 'slen', 'length','gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'score', 'pident', 'qseq', 'match','sseq']
	if overwrite or not os.path.isfile( out+'/result.xml' ):
		os.system( cmdstring )
	# parse the result to the filter
	result_handle = open('{o}/result.xml'.format(o=out))
	blast_records = NCBIXML.parse( result_handle )
	xml_list = []
	for rec in blast_records:
		for alignment in rec.alignments:
			for hsp in alignment.hsps:
				qacc = rec.query.split()[0]
				qacc = qacc.split('.')[0]
				if '|' in qacc:
					qacc = qacc.split('|')[1]
				sacc = alignment.accession.split()[0]
				sacc = sacc.split('.')[0]
				if '|' in sacc:
					sacc = sacc.split('|')[1]
				xml_list.append([qacc, sacc, rec.query_length, alignment.length, hsp.align_length, hsp.gaps, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, hsp.expect, hsp.score, 100*float(hsp.identities)/float(hsp.align_length), hsp.query, hsp.match, hsp.sbjct])
	df = pd.DataFrame(xml_list, columns=col_list)
	# calculate the query/subject coverage and add column
	'''df['sacc'] = df['sacc'].apply(lambda x: x.split()[0])
	df['qacc'] = df['qacc'].apply(lambda x: x.split()[0])
	df['sacc'] = df['sacc'].apply(lambda x: x.split('.')[0])
	df['qacc'] = df['qacc'].apply(lambda x: x.split('.')[0])'''
	df['pident'] = df['pident'].apply(lambda x: int(x))
	qcovs = []
	scovs = []
	for index, row in df.iterrows():
		# calculate and insert corresponding values
		qc = ((row["qend"]-row["qstart"]+1)*100)/row["qlen"]
		qcovs.append(qc)
		sc = ((row["send"]-row["sstart"]+1)*100)/row["slen"]
		scovs.append(sc)
	df.insert(loc=12, column='qcovs', value=qcovs)
	df.insert(loc=13, column='scovs', value=scovs)
	blast_filter( df, coverage, which, e_filter, alignment_length )
	# sort entries based on e-value, coverages
	df.sort_values(["evalue","qcovs"], inplace = True, ascending=[True,False])

	# try to run pfam for these results; extract sequence of good hits into one single fasta file
	genome_list = list(set(df['qacc'].tolist()))
	NCBI_list = list(set(df['sacc'].tolist()))
	# check the number of protein in the list, if empty, exit
	if len(genome_list) == 0:
		print >> sys.stderr, "No good results! Try to use other keywords or lower the filter standard."
		sys.exit(1)

	list_to_file( '{o}/good_genome.txt'.format( o=out ), genome_list )
	list_to_file( '{o}/good_NCBI.txt'.format( o=out ), NCBI_list )
	cmdstring = 'blastdbcmd -db {o}/blastdb/genome_blastdb -entry_batch {o}/good_genome.txt -out {o}/good_result.fasta;blastdbcmd -db {o}/blastdb/NCBI_blastdb -entry_batch {o}/good_NCBI.txt >> {o}/good_result.fasta'.format( o=out )
	os.system( cmdstring )
	# run hmmscan and parse the result as a table
	if overwrite or not os.path.isfile( out+'/pfam.out' ):
		cmdstring = 'hmmscan --cpu 4 --noali --cut_ga -o /dev/null --domtblout {o}/pfam.out /ResearchData/pfam/pfamdb/Pfam-A.hmm {o}/good_result.fasta'.format( o=out )
		os.system( cmdstring )
	# parse the result in pandas
	hmm_object = hmmParser( '{o}/pfam.out'.format( o=out ) )
	hmm_object.filterByCoverage(dcov)
	df_pfam = pd.DataFrame( hmm_object.matrix )
	domtblout_cols = 'target_name t_accession tlen query_name accession qlen evalue socre bias # of cevalue ievalue score bias hmm_from hmm_to ali_from ali_to env_from env_to acc description_of_target'.strip().split(' ')
	df_pfam.columns = domtblout_cols
	# get the clan accessions and clan info for the existing results
	df_pfam.insert( loc=3, column='clan_acc', value='' )
	df_pfam.insert( loc=4, column='clan_info', value='' )
	df_pfam['query_name'] = df_pfam['query_name'].apply(lambda x: x.split('.')[0])
	df_pfam['t_accession'] = df_pfam['t_accession'].apply(lambda x: x.split('.')[0])

	for index, row in df_pfam.iterrows():
		cmd = 'zgrep {s} /ResearchData/pfam/download/Pfam-A.clans.tsv.gz'.format( s=row["t_accession"] )
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
		output = process.stdout.read().split()
		#print output
		if 'CL' in output[1]:
			row["clan_acc"] = output[1]
			row["clan_info"] = output[2]
		else:
			row["clan_acc"] = 'N/A'
			row["clan_info"] = 'N/A'
	# for each pair of columns in the filtered output, map the value of each cell to the pfam result, check if the domians are the same
	domain_dic = clan_to_dic( 'query_name', 'target_name', df_pfam )
	clan_dic = clan_to_dic( 'query_name','clan_acc', df_pfam )
	
	#pprint.pprint(clan_dic)
	# create a list for unqualified results and drop rows according to the list in two dataframes
	bad_hits_blast = []
	#bad_hits_pfam = []
	no_domain = []
	for index, row in df.iterrows():
		'''if '|' in row["qacc"]:
			row["qacc"] = row["qacc"].split('|')[1]
		if '|' in row["sacc"]:
			row["sacc"] = row["sacc"].split('|')[1]'''
		# use set to compare if 2 lists have at least one common element
		try:
			if not set(domain_dic[row["qacc"]]) & set(domain_dic[row["sacc"]]):
				if not set(clan_dic[row["qacc"]]) & set(clan_dic[row["sacc"]]):
					df = df.drop(index)
					#bad_hits_pfam.append( row["qacc"] )
					#bad_hits_pfam.append( row["sacc"] )
		except KeyError:# in the case that the protein has no protein domains
			if not domain_dic.get(row["qacc"]):
				df_pfam.loc[-1]=['N/A', 'N/A', 'N/A', 'N/A', 'N/A', row["qacc"], 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A','N/A','N/A','N/A']
				df_pfam.index = df_pfam.index+1
				if row["qacc"] not in no_domain:
					no_domain.append(row["qacc"])
			if not domain_dic.get(row["sacc"]):
				df_pfam.loc[-1]=['N/A', 'N/A', 'N/A', 'N/A', 'N/A', row['sacc'], 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A','N/A','N/A','N/A']
				df_pfam.index = df_pfam.index+1
				if row["sacc"] not in no_domain:
					no_domain.append(row["sacc"])

	# delete rows according two lists
	df_pfam = df_pfam.drop_duplicates()
	df_pfam = df_pfam.set_index( 'query_name', drop=True)
	#df_pfam = df_pfam.drop(bad_hits_pfam) # drop by rows

	# save new sepeate results
	df['evalue'] = df['evalue'].apply(lambda x: '%.2e'%x)
	df.to_csv('{o}/good_blast_results.tsv'.format(o=out), sep='\t')
	df_pfam.to_csv('{o}/good_pfam_results.tsv'.format(o=out), sep='\t')
	pfam_dic = pfam_to_dic( df_pfam )
	# create hydropathy plots for these results, and add protein domains
	plot_path = out+'/plots'
	os.system( 'mkdir {p}'.format( p=plot_path ) )
	# create an final result and in the same loop finish data insertion
	outfile = open( out+'/result.html', 'w' )
	outfile.write( '<!DOCTYPE html><html><head><title>Missing Components Results</title><style type="text/css"></head><body>\n.label {text-align: right;width:50px;}\n.data {text-align:left;padding-left: 8px;width:100px;}\n.seq{border:2px solid black;height:70px; width:100%;overflow-x:auto;overflow-y:auto;margin:1em 0;background:grey;color:white;}tab1 { padding-left: 4em;}</style></head><body>\n' )
	# for each pair in the blast results (good one) create two seperated hydropathy plots and create a combined one with pfam domain covered, use blastdbcmd to extract sequences from the database and pass to the quod.py
	row_count = df.shape[0]
	if row_count < display:
		display = row_count
	df = df.head(display)
	df.reset_index(inplace = True)
	df.to_csv('{o}/good_blast_results1.tsv'.format(o=out), sep='\t') 
	for index, row in df.iterrows():
		d_str = domain_string( row["qacc"],row["sacc"], pfam_dic, no_domain )

		seq_string1 = 'blastdbcmd -db {o}/blastdb/genome_blastdb -entry {s}'.format( o=out, s=row['qacc'] )
		seq_string2 = 'blastdbcmd -db {o}/blastdb/NCBI_blastdb -entry {q}'.format( o=out, q=row['sacc'] )
		process = subprocess.Popen( seq_string1, stdout=subprocess.PIPE, shell=True )
		seq1 = process.stdout.read()
		process = subprocess.Popen( seq_string2, stdout=subprocess.PIPE, shell=True )
		seq2 = process.stdout.read()
		# create aligner
		alignment = row['qseq']+'\n'+row['match']+'\n'+row['sseq']

		# draw seperated plots first, draw bars of alignment and commnon domain
		# get the seperate command string first. It should have multiple domains and for same domain they should contain the same color
		os.system( 'quod.py -l {q} -q -s "{s1}" --width 15 -c red -d {p} -o {q}_vs_{s}.png -w {qs}-{qe} {ds}'.format( s1=seq1, p=plot_path, q=row['qacc'], qs=row['qstart'], qe=row['qend'], ds=d_str[0], s=row['sacc'] ) )
		os.system( 'quod.py -l {s} -q -s "{s2}" --width 15 -c blue -d {p} -o {s}_vs_{q}.png -w {ss}-{se} {ds}'.format( s2=seq2, p=plot_path, s=row['sacc'], ss=row['sstart'], se=row['send'], ds=d_str[1], q=row['qacc'] ) )
		# draw the aligned part as the combined graph
		seq = '>\n'+row['qseq']+'\n>\n'+row['sseq']
		temp = tempfile.NamedTemporaryFile(delete=True)
		try:
			temp.write(seq)
			temp.flush()
			os.system( 'quod.py {f} -l {q} -q --width 15 -d {p} -o {q}.png'.format( f=temp.name, p=plot_path, q=row['qacc']+'_'+row['sacc']+'_aligned' ) )
		finally:
			temp.close()
		
		# insert the blast info
		outfile.write( '<br /><hr style="border-style:solid; border-width:5px; color:black;"/><h2 style="text-align:left;">{q}</h2></font><font size = "4"><b>Hit Accession:</b>{s}</font><table width="600px" border="0" cellspacing="0" cellpadding="2"><tr><td class="label"><b>E-value:</b></td><td class="data">{evalue}</td><td class="label"><b>Identity:</b></td><td class="data">{pident}%</td><td class="label"><b>Length:</b></td><td class="data">{length}</td></tr><tr><td class="label"><b>Q_cov:</b></td><td class="data">{qcov}%</td><td class="label"><b>S_cov:</b></td><td class="data">{scov}%</td><td class="label"></td><td class="data"></td></tr></table><p><b>Alignment:</b><tab1>Query:{qstart}-{qend}<tab1>Subject:{sstart}-{send}</p><div class="seq"><pre>{align}</pre></div>'.format( q=row['qacc'],s=row['sacc'],length=row['length'],evalue=row['evalue'],pident=row['pident'],align=alignment, qcov=row['qcovs'],scov=row['scovs'], qstart=row['qstart'],qend=row['qend'],sstart=row['sstart'],send=row['send'] ) )
		# insert images to the html file
		outfile.write( '<center><table style = "width:100%" border = "0"><tr><td><center><a href = "{p}/{q}.png" target="_blank"><img src = "{p}/{q}.png" style="width:90%; height:90%"></a></center></td><td><center><a href = "{p}/{s}.png" target = "_blank"><img src= "{p}/{s}.png" style="width:90%; height:90%"></a></center></td></tr><tr><td colspan = "2"><center><a href = "{p}/{qs}.png" target = "_blank"><img src = "{p}/{qs}.png" style="width:50%; height:50%"></a></center></td></tr></table></center>'.format( p=os.path.join('./plots'), q=row['qacc']+'_vs_'+row['sacc'], s=row['sacc']+'_vs_'+row['qacc'],qs=row['qacc']+'_'+row['sacc']+'_aligned' ) )
		# insert the hmm info
		outfile.write( '<center><table style = "width:100%" border = "1"><tr><td>Domain</td><td>Domain_acc</td><td>Domain_len</td><td>Protein_acc</td><td>Protein_len</td><td>evalue</td><td>from</td><td>to</td><td>Clan</td><td>Clan_acc</td></tr>' )
		for obj in pfam_dic[row["qacc"]]:
			outfile.write( '<tr><td>{domain}</td><td>{dacc}</td><td>{dlen}</td><td>{pacc}</td><td>{plen}</td><td>{evalue}</td><td>{f}</td><td>{t}</td><td>{clan}</td><td>{cacc}</td></tr>'.format(domain=obj[6],dacc=obj[2],dlen=obj[7],pacc=obj[8],plen=obj[9],evalue=obj[10],f=obj[0],t=obj[1],clan=obj[11],cacc=obj[3]) )
		for obj in pfam_dic[row["sacc"]]:
			outfile.write( '<tr><td>{domain}</td><td>{dacc}</td><td>{dlen}</td><td>{pacc}</td><td>{plen}</td><td>{evalue}</td><td>{f}</td><td>{t}</td><td>{clan}</td><td>{cacc}</td></tr>'.format(domain=obj[6],dacc=obj[2],dlen=obj[7],pacc=obj[8],plen=obj[9],evalue=obj[10],f=obj[0],t=obj[1],clan=obj[11],cacc=obj[3]) )
		outfile.write( '</table></center><br>' )

	# Eventually create the result.html file
	outfile.write( '</body></html>' )
	outfile.close()


'''
 * Function name: domain_string( qlist, slist, df )
 * Description: this function is specific for formatting a cmdstring of quod.py. It passes lists of domains of both the query and the subject. for each common domain they share, it should be plotted in the same depth with same color.
 * Parameter: query list and subject list, which are keys and values in the dictionary generated previously, and the pfam dataframe
 * Return type: a list of cmdstrings containing both query string and subject string
'''
def domain_string( qacc, sacc, pfam_dic, except_list ):
	result = []
	#print except_list
	# special case: both qacc and sacc have no domains at all
	if any(qacc in s for s in except_list) and any(sacc in s for s in except_list):
		result.append('')
		result.append('')
		return result
	elif any(qacc in s for s in except_list):
		result.append('')
		string = Domain_string().one_protein(sacc,pfam_dic)
		result.append(string)
		return result

	elif any(sacc in s for s in except_list):
		string = Domain_string().one_protein(qacc,pfam_dic)
		result.append(string)
		result.append('')
		return result
	
	else:
		#print qacc
		#print sacc
		return Domain_string().two_proteins(qacc,sacc,pfam_dic)

'''
 * Function name: filter( blast_df, coverage, which, evalue )
 * Description: filter the blast result by evalue and coverage( query, subject, both, either )
 * Parameter: some parameters of get_ouput
 * return type: none
'''
def blast_filter( df, coverage, which, evalue, alignment_length ):
	#filter according to default setting or user specification
	df.ix[~(df['evalue'].astype(float) > evalue)]
	if which == 'b':
		df.ix[~(df['qcovs'] < coverage)]
		df.ix[~(df['scovs'] < coverage)]
	elif which == 'q':
		df.ix[~(df['qcovs'] < coverage)]
	elif which == 's':
		df.ix[~(df['scovs'] < coverage)]
	else:
		if df['qlen'].mean() < df['slen'].mean():
			df.ix[~(df['qcovs'] < coverage)]
		else:
			df.ix[~(df['scovs'] < coverage)]

	df.ix[~((df['length']-df['gaps']) < alignment_length)] 

	#df[(df['pident'] >= default_cover) & (float(df['evalue']) <= default_filter)] # the result shall only contain desired raws
'''
 * Function name: table_to_dic( col1, col2, pfam_df )
 * Description: this fucntion will convert a table to a dictionary containing a dictionary. The first layer is protein accessions; the second layer is domains that each protein has; the last layer is coordinates of each domain
 * Parameter: the name of the table file
 * Return type: the dictionary object. This function is specific to pfam dataframe
'''
def pfam_to_dic( df ):
	result = {}
	for index, row in df.iterrows():
		if index not in result:
			result[index] = [[row["env_from"],row["env_to"],row["t_accession"],row["clan_acc"],0,0, row["target_name"],row["tlen"],index,row["qlen"],row["evalue"],row["clan_info"]]]
		else:
			result[index].append( [row["env_from"],row["env_to"],row["t_accession"],row["clan_acc"],0,0, row["target_name"],row["tlen"],index,row["qlen"],row["evalue"],row["clan_info"]] )
	return result

def clan_to_dic( col1, col2, df ):
	result = {}
	for index, row in df.iterrows():
		if row[col1] not in result:
			result[row[col1]] = [row[col2]]
		else:
			result[row[col1]].append(row[col2])
	return result

'''
 * Function name: list_to_file(name, myList=[], *args)
 * Description: this function convert a list to a txt file with user-specified name
 * parameter: name indicating the filename
 * Return type: void
'''
def list_to_file( name, myList=[], *args ):
	data = open( name, 'w' )
	for accession in myList:
		data.write('{a}\n'.format( a=accession ))
	data.close()
	 

def main():
	# parse arguments
	parser = argparse.ArgumentParser()
	# create a parser to handle options and arguments. Show help info if no args
	parser.add_argument( '-k', '--keyword', type = str, dest = 'keywords',required = True, metavar = '<keywords string>', help = 'MANDATORY. Descriptions of the missing components for searching on the NCBI website. Logical expressions like "AND","OR" can be used for detailed description.' )
	parser.add_argument( '-l', '--local', action = 'store_true', dest = 'local', default = False, help = 'allow user to choose if the sequences are retrieved from a local database. The default setting is to download from NCBI online. If the number of sequences is much fewer than what were found online, it will retrieve from the local database. The final result will be the one with more sequences.' )
	parser.add_argument( '-g', '--genome', type = str, dest = 'genome',required = True, metavar = '<string input_file>', help = 'MANDATORY. Name or location of a complete genome file in FASTA format or FASTA compressed format(bz,gz).' )
	parser.add_argument( '-e', '--evalue', type = float, dest = 'evalue', metavar = '<evalue>', default = 1e-6, help = 'the e-value threashold for blast. Default value is 1e-6. Users can lower the e-value if they found nothing from blast.' )
	parser.add_argument( '-cf', '--coverage_filter', type = float, dest = 'coverage', metavar = '<a number between 0-100>', default = 50, help = 'percent coverage of the alignment of the smaller protein, used for filtering blast results. default is 50' )
	parser.add_argument( '-w', '--which', type = str, dest = 'which', metavar = '<x,q,s,b>', default = 'x', help = 'decide whether coverage filter will be applied on query(q), subject(s), both(b) or either(x) one. default is "x"') 
	parser.add_argument( '-ef', '--evalue_filter', type = float, dest = 'evalue_filter', metavar = '<evalue threashold for the result>', default = 0.000001, help = 'an additional filter for blast results. This filter ensures that e-values of all blast results must be equal to or smaller than the threshold. Default threshold is 1e-6' )
	parser.add_argument( '-dc', '--domain_coverage', type = float, dest = 'dcov', metavar = '<a number between 0-100>', default = 35, help = 'a domain filter for pfam results. Default is 35 percent. Domains which coverages are over the threshold will be considered as true hits.')
	parser.add_argument( '-al', '--alignment_length', type = int, dest = 'alignment_length', metavar = '<minimum alignment kept by the blast result>', default = 30, help = 'a filter for blast results that keeps alignment length (without gaps) greater than the threshold. Default value is 30')
	parser.add_argument( '-r', type = float, dest = 'redundancy', metavar = '<sequence identity threashold>', default = 0.9, help = 'a global sequence identity threashold for subject sequences. Default is 0.9. It removes redundant sequences whose similarities are higher than 0.9' )
	parser.add_argument( '-d', '--display', type = int, dest = 'display', default = 100, help = 'the minimum number of results displayed in the final html output. Default is 100. if the number of blast result is less than the threshold, it will display all results')
	parser.add_argument( '-out', type = str, dest = 'out', metavar = '<string output_directory>', default = 'MissingComponent_result', help = 'the name, or path and name of the output directory containing all files generated by this program and the final result in html format. Default name is "MissingComponent_result" generated in the current directory' )
	parser.add_argument( '-o', '--overwrite', action = 'store_true', dest = 'overwrite', default = False, help = 'overwrite all output files. Disabled by default' )
	args = parser.parse_args()
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	# do actions about the overwrite option
	args.dcov = args.dcov/100
	location = os.getcwd() + '/' + args.out
	if not args.overwrite:
		if os.path.isfile( location+'/result.html' ):
			sys.exit()
	# check if the genome file is valid or not
	if not os.path.isfile( args.genome ):
		print >> sys.stderr, "Genome file does not exist!"
		sys.exit(1)
	else:
		if ".bz2" in args.genome:
			genomeFile = 'bzip2 -c {f} | '.format( f=args.genome )
		elif ".gz" in args.genome:
			genomeFile = 'gunzip -c {f} | '.format( f=args.genome )
		else:
			genomeFile = 'cat {f} | '.format( f=args.genome )

	# retrieve sequences from NCBI and create the output directory, create the blast databse at the same time
	if retrieve_seq( args.keywords, args.out, args.local, args.overwrite, args.redundancy ):
		dbname = create_db( genomeFile, args.out )
	else:
		print >> sys.stderr, "No results from NCBI! Try to use other keywords."
		sys.exit(1)

	if dbname == 'NCBI':
		get_output( args.display, args.dcov, args.overwrite,genomeFile, location+'/blastdb/NCBI_blastdb', args.evalue, args.coverage, args.which, args.evalue_filter, location, args.alignment_length)
	else:
		get_output( args.display, args.dcov, args.overwrite,location+'/NCBI_seq.fasta', location+'/blastdb/genome_blastdb', args.evalue, args.coverage, args.which, args.evalue_filter, location, args.alignment_length)


main()
