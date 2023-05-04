#!/usr/bin/env python

#imports parser,functions of the TMS Split program, and system module
#changes implemented in tmssplit function
import argparse
import tmsFunction as tms
import sys


#------------------ PARSER-----------------------#

def Main():

	parser = argparse.ArgumentParser(description="TMS Split is a program that allows for the cutting of sequences in FASTA format. The program has four operations: Residue Range Cut,Equal Cut, TMS Split, and TMS Range Cut.\n\n This program runs in Python3. \n\n  Example usage: tmsSplit -if /path/to/input/file -od /path/to/output/directory -non -tms 4 -append")

	parser.add_argument('-if',help='The path to the file containing the input sequence(s) in FASTA format.',dest='inputFile', type=str)
	parser.add_argument('-od',help='The path to the output directory.' +' If the directory is not present,it will be created.' , dest='outputDir', type=str)
	parser.add_argument('-of', help='Name of the output file in which the processed sequence will be placed. NOTE: Do not add any extensions to the name. ".faa" will automatically be assigned.',default = '', dest='outFile', type=str)

	operation = parser.add_mutually_exclusive_group()
	operation.add_argument('-rangeCut',action='store_const', const='0', dest='operation', help='Extract a segment of variable size from a protein sequence. Start(-s or --start) and End(-e or --end) required for range cut.')
	operation.add_argument('-equal', action='store_const', const='1',dest='operation',help='Cut protein(s) into segments of equal length. Parts (-p or --parts) required for equal cut.')
	operation.add_argument('-tmsCut', action='store_const',const='2',dest='operation', help="Cut multiple sequences between any two TMS's. Start(-s or --start) and End(-e or --end) required for TMS cut.The tail argument is optional, but the default is 3. ")
	operation.add_argument('-split', action='store_const', const = '3',dest='operation',help='Split protein(s) into groups of TMS. Choose overlapping(-over) or non-overlapping(-non). DEFAULT is non-overlapping. TMS(-tms) argument required for TMS split.The tail argument is optional, but the default is 3.')

	parser.add_argument('-t','--tail',type=int,dest='tail', default=3,help='The number of residues left before and after the sequence segments. DEFAULT is 3.')

	parser.add_argument('-p','--parts',dest='parts',type=int, help='Argument for equal cut. The number of segments of equal length the sequence(s) should be cut into.')

	parser.add_argument('-s','--start',dest='start',type=int, help='Argument for TMS cut. The number of the TMS or Residue that will be the first in the segment. ')
	parser.add_argument('-e','--end',dest='end',type=int, help='Argument for TMS cut. The number of the TMS or Residue that will be the last in the segment.')

	parser.add_argument('-tms',dest='tms',type=int, help='Argument for TMS Split. The number of TMS per group.')

	split = parser.add_mutually_exclusive_group()
	split.add_argument('-over', action='store_const',dest='overlap', const=1, help='Argument for TMS Split. TMS are grouped by overlapping sections. Example: [1,2,3],[2,3,4] etc...')
	split.add_argument('-non',action='store_const',dest='overlap',const=0, help='DEFAULT Argument for TMS Split. TMSs are grouped by non-overlapping sections. Example: [1,2,3],[4,5,6] etc... Also requires a condition -ignore, -append, or -new')

	remain = parser.add_mutually_exclusive_group()
	remain.add_argument('-ignore', action='store_const',const=0,dest='remainingSeq', help='Argument for TMS Split Non-Overlapping. Any TMS that are not able to be formed into a full group are ignored.')
	remain.add_argument('-append', action='store_const',const=1,dest='remainingSeq', help='DEFAULT Argument for TMS Split Non-Overlapping. Any TMS that are not able to be formed into a full group are added to the last  possible segment.')
	remain.add_argument('-new', action='store_const',const=2,dest='remainingSeq',help='Argument for TMS Split Non-Overlapping. Any TMS that are not able to be formed into a full group are added to another segment')


	#Retrieve Parser Values
	args = parser.parse_args()


	#Evaluate if input file and output directory are given
	if args.inputFile and args.outputDir:

#----------------------Residue Range Cut-------------------------#
		if args.operation == '0':

			# Test if the start and end arguments are provided and are greater than 0
			if args.start and args.end and args.start > 0 and args.end > 0:

				#assigning argument values to variables
				start = args.start
				end = args.end

				#Check whether start is less than end
				#Throw an error to user
				if start >= end:
					parser.error('Start value was greater than End value.')

				#execute functions for Residue Range Cut
				#Catch any Exceptions
				try:
					tms.inputRangeCut(args.inputFile,args.outputDir,args.outFile,start,end)
				except Exception as e:
					print(str(e))
			else:
				parser.error('-s and -e required for Cut Residue Range. Values must be greater than 0 and cannot be equal to each other.')


#---------------------Equal Cut------------------------------#
		elif args.operation == '1':

			#Check if the parts argument is selected and has a value greater than 0
			if args.parts and args.parts > 0:

				#Execute Equal Cut Operation
				try:
					tms.inputEqualSplit(args.inputFile,args.outputDir,args.outFile,args.parts)
				except Exception as e:
					print(str(e))

			else:
				parser.error('-p (integer > 0) reqiuired for -equal operation.')



#--------------------TMS Range Cut-------------------------#
		elif args.operation == '2':

			# Test if the start and end arguments are provided and are greater than 0
			if args.start and args.end and args.start > 0 and args.end > 0:

				#assigning argument values to variables
				start = args.start
				end = args.end

				#Check whether start is less than end
				#Throw an error to user.
				if start > end:
					parser.error('Start value was greater than End value.')

				#Execute TMS Cut operation
				try:
					tms.tmsCutInput(args.inputFile,args.outputDir,args.outFile,start,end,args.tail)
				except Exception as e:
					print(str(e))
			else:
				parser.error('-s and -e (integers > 0) required for -cut operation.')



#-------------------TMS Split--------------------------#
		elif args.operation == '3':

		#++++++++++++DEFAULTS++++++++++++#
			#overlapping and nonoverlapping(Non-overlapping is default)
			overlap = 0
			#Remaining Sequence if Condition is non-overlapping(append is default)
			remainingSeq = 1

		#++++++++If options are selected++++++++++#

			#Overlapping Condition:
			#Non-Overlapping = 0(default)
			#Overlapping = 1
			
			#Default Case (python3 change)
			if (args.overlap == None):
				args.overlap = 0
			if args.overlap >= 0:


				overlap = args.overlap


			#Remaining TMS for Non-Overlapping
			#ignore = 0
			#append = 1(default)
			#new = 2
			
			#Default Case (python3 change)
			if (args.remainingSeq == None):
				args.remainingSeq = 1
			if args.remainingSeq >= 0:
				
				remainingSeq = args.remainingSeq



			#Execute tms Split operation
			if args.tms and args.tms > 0:
				try:
					tms.tmsSplit(args.inputFile,args.outputDir,args.outFile,args.tms,args.tail,overlap,remainingSeq)
				except Exception as e:
					print(str(e))
			else:
				parser.error('TMS argument (-tms) required for split operation.It must be a positive integer')

	#No Operation Argument is given
	#So print help statement and exit
	else:

		parser.print_help()
		sys.exit(1)


Main()
