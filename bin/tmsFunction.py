# /Users/pranaviddamsetty/anaconda2/bin/python
#xrange changed to range

from Bio import SeqIO
import os.path
import random
import string
import re
import subprocess
from subprocess import *
import shlex

'''

'''

def inputRangeCut(file_path,outdir,outFile,start,end):

	checkInput(file_path,outdir)

	in_file = open(os.path.relpath(file_path),'rU')

	errors = ''
	for record in SeqIO.parse(in_file,"fasta"):
		sequence = record.seq
		seqID = record.id

		output =''
		if outFile != '':
			output = outFile
		else:
			seqName = seqID.replace(':','_')
			output = re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1)

		errors += RangeCut(sequence,seqID,start,end,outdir,output)

	raise Exception(errors)





def RangeCut(seq,seqID,start,end,outdir,outFile):

	if start < len(seq) and end <= len(seq):
		sequence = list(seq)

		fin_seq = ''

		fin_seq += '>{} Residue {} to {}\n'.format(seqID,start,end)

		fin_seq += '{}\n'.format("".join(sequence[start - 1: end ]))

		writeFile(fin_seq,outdir,outFile,'rangeCut')

		return ''
	else:

		seqName = seqID.replace(':','_')
		ID = re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1)

		return 'Sequence {} only has {} residues.\n'.format(ID, len(seq))


"""
Function Input recieves input from the user and evaluates parameters for validity.
Once all validations are met, parameters are sent to equalSplitSeq with individual
fasta sequences.
"""
def inputEqualSplit(file_path,outdir,outFile,parts):


	checkInput(file_path,outdir)

	in_file = open(os.path.relpath(file_path), "rU")

	#Prepare individual sequence and send to equalSplitSeq
	for record in SeqIO.parse(in_file, "fasta"):
		sequence = record.seq
		seqID = record.id

		output =''
		if outFile != '':
			output = outFile
		else:
			seqName = seqID.replace(':','_')
			output = re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1)

		equalSplitSeq(sequence,seqID,outdir,output,parts)

"""
checkInput validates existence of the input file. If the file does
not exist, then it returns an error message. If file does
exist, it returns a reference to the opened file.

It also validates the existence of the output directory.
If it does not exist, the output directory is created.
"""
def checkInput(file_path,outdir):

	#Validate existence of Input File
	if os.path.exists(file_path):
		pass
	else:
		raise Exception('Error: {} not found'.format(file_path))
		quit()

	#Validate existence of output directory(else create directory)
	if os.path.exists(outdir):
		pass
	else:
		#print('Directory Not Found: Creating Directory')
		os.makedirs(outdir)





"""
 Function equalSplitSeq takes a sequence and cuts it into the number of
 parts requested by the user. If the function cannot be cut equally into
 the number of parts requested, the remainder will be added to the final
 segment.
"""

def equalSplitSeq(str_seq,seqID,outdir,outFile,parts):
	#Prepare sequence parameters
	sequence = list(str_seq)
	length = len(sequence)
	segLength= int(length/parts)
	fin_seq = ''
	
	#Create segments and add to string 'fin_seq'
	for x in range(0,parts):
		if x==parts-1:
			end = length
			out_seq = ''.join(sequence[((parts-1)*segLength): ])
		else:
			end = (x*segLength)+segLength
			out_seq = ''.join(sequence[(x*segLength):end])

		fin_seq += ('>{} Residue {} to {} | Length: {}\n'.format(seqID,(x*segLength)+1,end,len(out_seq)))
		fin_seq += ("{}\n".format(out_seq))
	writeFile(fin_seq,outdir,outFile,'ecut_seq')
	return





"""
Function tmsSplit
	recieves:
	-file path( for input)
	-output directory
	-the number of tms per segments
	-the size of the tail(both before the segment and after)
	-the overlap attribute
		- non-overlapping
		- overlapping
	- the instructions for the remainng tms(only applies to non-overlapping)
		- ignore
		- append to the last segment
		- create a new segment
"""
def tmsSplit(file_path, outdir,outFile,tmsPerSeq,tail_size,overlap,remainSeq):

	checkInput(file_path, outdir)
    #opens file so that it is ready to read

	in_file = open(os.path.relpath(file_path), "r")

    #creates the shell command to execute HMMTOP

	cmnds = 'hmmtop -if={} -sf=FAS -pi=spred -is=pseudo 2> /dev/null'.format(file_path)


    #prepares arguments and executes HMMTOP
    # output is captured into variable hmm_out

	args = shlex.split(cmnds)
	output = subprocess.Popen(args,shell=False, stdout=subprocess.PIPE, stderr = subprocess.PIPE,universal_newlines=True)
	result, error = output.communicate()
	hmm_out = result


    #Creates two lists: one containing the sequences and one containing their fasta formated headers

	sequences = []
	records = []

	for record in SeqIO.parse(in_file, "fasta"):
		sequences.append(record.seq)
		records.append(record.id)

	#reads result of hmmtop and determines number
	#and location of tms
	lines = str(hmm_out).rstrip().split('\n')
	tms_list = [[]]



    #cleans the output from HMMTOP and creates lists containing tms pairs
	for line in lines:
		if not bool(line) and bool(re.search('^\>HP.+',line)):
			continue
		tms_line = re.split('\s+(IN|OUT)\s+',line)
		tmsNumbers = tms_line[2].rstrip().split()
		if tmsNumbers[0] == 0:
			continue
		tms_index = hmmPair(tmsNumbers)

		tms_list.append(tms_index)


    #initializes a string that captures all error statements
	errors = ''

    # sends individual sequences to the relavent operation: either overlapping or non-overlapping
    # Also initializes file name for each sequence so that it may be identified
	for x in range(0, len(sequences)):
		sequence = sequences[x]
		record = records[x]
		tms_index = tms_list[x+1]
		output = ''
		if outFile != '':
			output = outFile + '_' + str(x+1)
		else:
			seqName = record.replace(':','_')
			output = re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1)
		if overlap == 0:
			errors += nonOverlapTMS(sequence,record,tms_index,tmsPerSeq,tail_size,remainSeq,outdir,output)
		if overlap == 1:
			errors += overlapTMS(sequence,record,tms_index,tmsPerSeq,tail_size,outdir,output)

    #collectively raises any errors that may have occured
	raise Exception(errors)
	



"""
    nonOverlapTMS is one of two helper functions for the TMS split operation
    It recieves:
    
        -an input sequence
        -the record id for the sequence
        -a list of the tms in the sequence
        - the number of tms per group
        - the tail size
        -the condition for the remaining sequences ( whatever TMS are not evenly divided)
            -ignore
            - append to the last formed group
            - append to a new group
        - the output direcotry
        -the  output file name
        
"""
def nonOverlapTMS(sequence,record,tms_index,tmsSeq,tail,remainSeq,outdir,outFile):

	seq = list(sequence)
	length = len(seq)
	numTMS = int(tms_index[0])
	segments = numTMS/tmsSeq
	remainTMS = numTMS % tmsSeq
	finalSeq =''

	if numTMS != 0:
		for x in range (0,int(segments)):
			first_tms = tms_index[(tmsSeq*x)+1]
			last_tms = tms_index[(tmsSeq*x)+tmsSeq]
			start = int(first_tms[0])
			startIndex = tmsTailStart(tail,start)

			end = int(last_tms[1])
			endIndex = tmsTailEnd(tail,length,end)

			#append to last segment
			if x == segments - 1 and remainSeq == 1:

				last_tms = tms_index[-1]
				end = int(last_tms[1])
				endIndex = tmsTailEnd(tail,length,end)

				finalSeq += ('>{} | TMS {} to {}\n'.format(record,(tmsSeq*x)+1,numTMS))

				finalSeq +=('{}\n'.format(''.join(seq[startIndex:endIndex])))
		
            #create a new segment
            #-------Bug Fixed: Now Checks if there are tms remaining------------
			elif x == segments - 1 and remainSeq == 2 and remainTMS >= 1:

				finalSeq += ('>{} | TMS {} to {}\n'.format(record,(tmsSeq*x)+1,(tmsSeq*x)    +tmsSeq))

				finalSeq +=('{}\n'.format(''.join(seq[startIndex:endIndex])))

				first_tms = tms_index[(tmsSeq*x)+tmsSeq+1]
				start_final = int(first_tms[0])
				startFinal = tmsTailStart(tail,start_final)

				last_tms = tms_index[-1]
				final = int(last_tms[1])
				endFinal = tmsTailEnd(tail,length,final)

				finalSeq += ('>{} | TMS {} to {}\n'.format(record,(tmsSeq*x)+tmsSeq+1,numTMS))
				finalSeq +=('{}\n'.format(''.join(seq[startFinal: endFinal])))


			#ignore
			else:

				finalSeq += ('>{} | TMS {} to {}\n'.format(record,(tmsSeq*x)+1,(tmsSeq*x)+tmsSeq))
				finalSeq +=('{}\n'.format(''.join(seq[startIndex:endIndex])))


		writeFile(finalSeq,outdir,outFile,'nosplit_tms')

		return ''

	else:

		seqName = record.replace(':','_')
		return 'Sequence {} has 0 TMS\n'.format(re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1))






"""
Function for overlapping condition for TMS Split.
"""
def overlapTMS(sequence,record,tms_index,tmsSeq,tail,outdir,outFile):
	seq = list(sequence)
	length = len(sequence)
	numTMS = int(tms_index[0])

	finalSeq = ''

	if numTMS != 0:

		for x in range(0,numTMS-tmsSeq+1):
			first_tms = tms_index[x+1]
			last_tms = tms_index[x+tmsSeq]

			start = int(first_tms[0])
			startIndex = tmsTailStart(tail,start)

			end = int(last_tms[1])
			endIndex = tmsTailEnd(tail,length,end)

			finalSeq += ('>{} | TMS {} to {}\n'.format(record,x+1,x+tmsSeq))
			finalSeq += ('{}\n'.format(''.join(seq[startIndex:endIndex])))

		writeFile(finalSeq,outdir,outFile,'osplit_tms')

		return ''

	else:

		seqName = record.replace(':','_')
		return 'Sequence {} has 0 TMS\n'.format(re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1))




"""
Recieves the hmmtop results and creates pair listing the
start and end of the tms
"""
def hmmPair(tms_info):

	numTMS = tms_info[0]

	tms = [[]]

	tms[0]=numTMS
	for x in range(0, int( len(tms_info)/2)):
		tempList = [ tms_info[(2*x)+1],tms_info[(2*x)+2]]
		tms.append(tempList)

	return tms



"""
Accepts input for the tmsCut operation
"""

def tmsCutInput(file_path, outdir,outFile,first,last,tail):

	checkInput(file_path,outdir)

	in_file = open(os.path.relpath(file_path),mode = "r",newline = None)


	sequences = []
	records = []

	for record in SeqIO.parse(in_file, "fasta"):
		sequences.append(record.seq)
		records.append(record.id)

	cmnds = 'hmmtop -if={} -sf=FAS -pi=spred -is=pseudo 2> /dev/null'.format(file_path)


	args = shlex.split(cmnds)
	output = subprocess.Popen(args,shell=False, stdout=subprocess.PIPE, stderr = subprocess.PIPE,universal_newlines=True)
	hmm_out, error = output.communicate()
	lines =hmm_out.rstrip().split('\n')
	tms_list = [[]]
	for line in lines:
		if not (bool(line) and bool(re.search('^\>HP.+',line))):
			continue
		tms_line = re.split('\s+(IN|OUT)\s+',line)
		tmsNumbers = tms_line[2].split()
		if tmsNumbers[0] == 0:
			continue
		tms_index = hmmPair(tmsNumbers)
		tms_list.append(tms_index)

	errors = ''
	for x in range(0, len(sequences)):
		sequence = sequences[x]
		record = records[x]
		tms_index = tms_list[x+1]

		output = ''
		if outFile != '':
			output = outFile + '_' + str(x+1)
		else:
			seqName = record.replace(':','_')
			output = re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1)

		errors += tmsCut(sequence,record,tms_index,first,last,tail,outdir,output)
	raise Exception(errors)



"""
tmsCut is the function which cuts the sequences from one tms to another
	recieves: sequence, record, tms locations,first tms, last tms, tail
"""

def tmsCut(sequence,record,tms_index,first,last,tail,outdir,outFile):

	seq = list(sequence)
	length = len(seq)
	numTMS = int(tms_index[0])

	finalSeq =''

	if first > 0 and last <= numTMS:

		first_tms = tms_index[first]
		last_tms = tms_index[last]

		start = int(first_tms[0])
		startIndex = tmsTailStart(tail,start)

		end = int(last_tms[1])
		endIndex = tmsTailEnd(tail,length,end)

		finalSeq += ('>{} | TMS {} to {}\n'.format(record,first,last))
		finalSeq += ('{}\n'.format(''.join(seq[startIndex:endIndex])))

		writeFile(finalSeq,outdir,outFile,'tms_range')

		return ''
	else:

		seqName = record.replace(':','_')
		return('{}: TMS out of Bound. Only {} TMS in sequence.\n'.format(re.search(r'([0-9a-zA-Z_\.-]+)',seqName).group(1),numTMS))



"""
    
    
"""
def tmsTailStart(tail,tms):

	if tail >= tms:

		index = tms - 1

		return index
	else:

		return tms-tail-1


"""
"""
def tmsTailEnd(tail,length,tms):

	if tail >= length-tms:

		index = length - 1

		return index

	else:

		return tms+tail-1




"""
writeFile creates the file and writes the final sequence
onto the file.
"""

def writeFile(sequence,outdir,outFile,operation):

	fileName = outdir + '/' + outFile + '.faa'

	out_file = open(fileName, 'w')
	out_file.write(sequence)
	out_file.close()
