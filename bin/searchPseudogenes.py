#!/usr/bin/env python

'''
 * File name: Blastx.py
 * Author: Yichi Zhang
 * Date: 2018/2/22
 * Description:  The program is a interface of blastx in the lab of Dr.Saier, written by Pyhton 2.7.
 * 		 It aims to automatically get protein sequence and run ORF scanning to get results by entering the transp 
 *		 -orter system accession and specific proteinidentifer.
 * Parameter: flags - let users to choose which output to form
 * 	      DNA sequence - string of the complete DNA sequence FASTA file
 * 	      system - the accession of a transportor system
 *            protein - the protein identifier
 * Return type: none
 * Side-effects: form a output file of blastx results
'''

import os
import sys
import getpass
import subprocess
import bz2
import gzip
import argparse

'''
 * Function name: createdb()
 * Desription: to generate a tcdb on the current user account
 * Parameter: N/A
 * Return type: void
'''
def createdb():
	#create db and blastdb directory
	os.system ( 'mkdir ~/db' )
	os.system ( 'extractFamily.pl -i all -o ~/db/blastdb -f blast' )
	os.system ( 'extractFamily.pl -i all -o ~/db/blastdb -f fasta' )

'''
 * Function Name: createDirectory( filename )
 * DEscription: check if filename exists as a directory in the current location. If not, create a new one
 * Parameter: string - filename
 * Return Type: void
'''
def createDirectory( filename ):
	if not os.path.isdir( filename ):
		cmdstring = 'mkdir {f}'.format( f = filename )
		os.system( cmdstring )

'''
 * Function Name: checkProtein( filename, directory )
 * Description: check if specific protein exists in the location directory.
 * Parameter: filename - the protein interested; directory - location directory
 * Return type: True - if the protein file does not exist; False - the protein file has already exists
'''
def checkProtein( filename, directory ):
	result = False
	location = './' + directory
	if not any( File.startswith(filename) for File in os.listdir( location )):
		result = True
	return result

'''
 * Function Name: getProtein( accessions, pname, directory )
 * Decription: this function calls blastdbcmd and extract protein sequences by accessions into an output file
 * Paramenters: a list of protein identifiers; a directory name for output location; the name of output protein file
 * Return type: string of location of the protein file
 * Side effect: generate a FASTA file fot that sequence
'''
def getProtein( accessions, pname, directory ):
	result = directory + '/' + pname + '.faa'
	cmdstring = 'blastdbcmd -db tcdb -entry {protein} -out {outfile}'.format( protein = accessions, outfile = result )
        os.system(cmdstring)

	return result
'''
 * Function Name: getProteinFromFile( filename, pname, directory )
 * Description: this function calls blastdbcmd and extract protein sequences by an accession file into an output file
 * Parameters: accession filename; protein name for output file; directory name for output file
 * Return type: a string of location of the protein file for blast
 * Side effect: generate a fasta file
'''
def getProteinFromFile( filename, pname, directory ):
	result = directory + '/' + pname + '.faa'
	cmdstring = 'blastdbcmd -db tcdb -entry_batch {protein} -out {outfile}'.format( protein = filename, outfile = result )
	os.system( cmdstring )

	return result
'''
 * Function Name: checkOutput( option, location )
 * Description: this function checks if the location directory contains the correct output.
 * Parameters: option - the correct form of output; location - the directory that contains the output
 * Return type: True - if the directory has no correct output
		False - if the directory has output that people want
 * Side effect - None
'''
def checkOutput( option, location ):
	location = './' + location
	result = False
	if option == 'table':
		if not any( File.startswith('table') for File in os.listdir( location )):
			result = True
	elif option == 'alignment':
		if not any( File.startswith('alignment') for File in os.listdir( location )):
			result = True
	elif option == 'both':
		if not any( File.startswith('table') for File in os.listdir( location )) and not any( File.startswith('alignment') for File in os.listdir( location )):
			result = True
	return result

'''
 * Function name: getOutput( query:str, subject:str, evalue:float, comp_based_stats:str, seg:str, out:str, outfmt:str )
 * Description: this function will finally call blastx and the shell and get the final output
 * Parameters: all parameters for the blast
 * Return type: void
 * Side effect: create the output file for blast results
'''
def getOutput( query, subject, evalue, comp_based_stats, low_comp_filter, out, outfmt):
	cbs = comp_based_stats
	lcf = low_comp_filter
	
	if comp_based_stats:
		cbs = '2'
	if low_comp_filter:
		lcf = 'yes'
	table = '{q}blastx -query - -out ./{o}/table.out -evalue {e} -comp_based_stats {c} -seg {s} -subject {sj} -use_sw_tback '.format( q=query, o=out, e=evalue, c=cbs, s=lcf, sj=subject)
	alignment = '{q}blastx -query - -out ./{o}/alignment.out -evalue {e} -comp_based_stats {c} -seg {s} -subject {sj} -use_sw_tback '.format( q=query, o=out, e=evalue, c=cbs, s=lcf, sj=subject)
	cmdstring2 = ''
	fmt1 = '-outfmt "7 sacc slen qframe qstart qend sstart send score evalue length pident"'
	fmt2 = '-outfmt "0"'
	
	if outfmt == 'table':
		cmdstring1 = table + fmt1
	elif outfmt == 'alignment':
		cmdstring1 = alignment + fmt2
	else:
		cmdstring1 = table + fmt1
		cmdstring2 = alignment + fmt2
	
	os.system( cmdstring1 + ';' + cmdstring2 )


def main():
	# parse the arguments
	parser = argparse.ArgumentParser()
	#create a parser to handle options and arguments. No args, show help
	parser.add_argument( '-g', '--genome', type = str, dest = 'genome', metavar = '<string input_file>', help = 'MANDATORY. Name or location of a complete DNA sequence file in FASTA format, bzip2 format, or gzip format' )
	parser.add_argument( '-p', '--protein', type = str, dest = 'query', metavar = '<query string>', help = 'MANDATORY. Provide one of following:  proteinsequence.faa - file containing protein sequence of interest in FASTA format.        Single or list of TCID-Accessions seperated by commas - eg: 1.A.1.1.1-P0A334,1.A.1.10.1-P08104.       Name of a file containing TCID-Accessions seperated by new lines.' )
	parser.add_argument( '-e', '--evalue', type = float, dest = 'evalue', metavar = '<evalue>', default = 1.0, help = 'the e-value threashold for blast. Default value is 1.0' )
	parser.add_argument( '-c', '--comp-based-stats', action = 'store_true', dest = 'cbs', default = '0', help = 'use composition-based statistics. Disabled by default' )
	parser.add_argument( '-l', '--low-comp-filter', action = 'store_true', dest = 'lcf', default = 'no', help = 'filter low complexity regions. Disabled by default' )
	parser.add_argument( '-out', type = str, dest = 'out', metavar = '<directory name>', default = 'pseudogenesResult', help = 'the name of output directory which contains all output files. Location with name of directory is also accepted. Default name is "pseudogenesResult" generated in the current directory' )
	parser.add_argument( '-f', '--format', type = str, dest = 'format', metavar = '<string format>', default = 'table', help = 'table - generate the table of results;alignment - generate paired alignments;both - generate both table and alignment output file.Default is table' )
	parser.add_argument( '-o', '--overwrite', action = 'store_true', dest = 'overwrite', default = False, help = 'overwrite all output files. Disabled by default' )
	args = parser.parse_args()
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	# do actions about the arguements 
        # check if the database exsists. No, generate one.
	username = getpass.getuser()
	dbPath = '/Users/{user}/db'.format(user=username)
	if not os.path.isdir(dbPath):
		createdb()
	
	# check the format of genome and do actions, extract compressed files or direct store
	if os.path.isfile( args.genome ):
		if ".bz2" in args.genome:
			genomeFile = 'bzip2 -c {f} | '.format( f=args.genome )
		elif ".gz" in args.genome:
			genomeFile = 'gunzip -c {f} | '.format( f=args.genome )
		else:
			genomeFile = 'cat {f} | '.format( f=args.genome )
	else:
		parser.print_help()
		sys.exit(1)

	# check if the subject is a file of protein sequence, an accession list, or an accession file; get the tranporter system identifier
	if os.path.isfile(args.query):
		# open the file and read lines to form a list of lines
		with open ( args.query ) as f:
			accessions = f.readlines()
		# check the content to decide if it is a fasta or accession list
		if '>' in accessions[0]:
			proteinFile = args.query
			familyid, dash, proteinid = args.query.partition('.') # get the transporter system id and protein id
			createDirectory( args.out ) # check and create a directory to hold all output if this directory does not exist
		else:
			accessions = [x.strip() for x in accessions]
			familyid, dash, proteinid = accessions[0].partition('-')
			createDirectory( args.out )
			# extract the protein sequence from database and put it into the directory created before
			if args.overwrite or checkProtein( proteinid, args.out ):
				proteinFile = getProteinFromFile( args.query, proteinid, args.out )
			else:
				proteinFile = args.out + '/' + proteinid + '.faa'
	else:
		# create a list for input query strings
		accessions = args.query.split(',')
		familyid, dash, proteinid = accessions[0].partition('-')
		createDirectory( args.out )
		if args.overwrite or checkProtein( proteinid, args.out ):
			proteinFile = getProtein( args.query, proteinid, args.out )
		else:
			proteinFile = args.out + '/' + proteinid + '.faa'

	# form the output
	# check if any output exists or overwirte is on 
	if args.overwrite or checkOutput(args.format, args.out) :
		getOutput( genomeFile, proteinFile, args.evalue, args.cbs, args.lcf, args.out, args.format )

main()
	
