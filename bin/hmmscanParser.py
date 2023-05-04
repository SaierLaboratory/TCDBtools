#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from re import split
from itertools import groupby
from operator import itemgetter

'''
 * file name: hmmParser.py
 * Description: this file suppose to parse the domain table of hmmscan and apply certain filters to the result. Special thanks to Enzo Guerrero-Araya for inspiration and coding help
'''
__author__ = 'Enzo Guerrero-Araya (biologoenzo@gmail.com), Yichi Zhang (yiz370@ucsd.edu)'
class hmmParser( object ):

	def __init__(self, hmmfile):
		self.hmmfile = hmmfile
		self.data_hmmfile = open(self.hmmfile).read()
		self.parameters = {}  # dict_keys(["Target file", "Option settings", "Program", "Version", "Date", "Current dir", "Pipeline mode", "Query file"])
		for i, line in enumerate(self.data_hmmfile.split("\n")[-10:-2]):
			line = line.split(":")
			key = line[0][2:]
			value = line[1].strip()
			self.parameters[key] = value
		self.matrix = self.hmmscanParser()


	def filterByEvalue(self, evalue=1e-18):
		for i, row in enumerate(self.matrix):
			if float(row[12]) > evalue:  # domain[11] -> i-Evalue (whole database)
				self.matrix.pop(i)

	def filterByBitscore(self, bits=50):
		for i, row in enumerate(self.matrix):
			if float(row[13]) < bits:  # domain[13] -> Bitscore
				self.matrix.pop(i)

	def filterByCoverage(self, cov=0.35):  # covered fraction of HMM
		if self.parameters["Program"] == "hmmscan":
			for i, row in enumerate(self.matrix):
				coverage = (float(row[16]) - float(row[15])) / float(row[2])
				if coverage < cov:  # domain[13] -> Bitscore
					self.matrix.pop(i)

	def uniqueByBestBitscore(self,):  # by domain
		matrix = sorted(self.matrix, key=itemgetter(3), reverse=True)
		for query, group in groupby(matrix, itemgetter(3)):
			group = list(group)
			if len(group) == 1:
				continue
			group = sorted(group, key=itemgetter(13), reverse=True)
			for dom in group[1:]:
				index = self.matrix.index(dom)
				self.matrix.pop(index)

	def hmmscanParser(self, ):
		matrix = []
		# header = ["target name", "accession", "tlen", "query name",
		#           "accession", "qlen", "E-value", "score", "bias", "#", "of",
		#           "c-Evalue", "i-Evalue", "score", "bias", "from", "to",
		#           "from", "to", "from", "to", "acc", "description of target"]
		# matrix.append(header)
		for line in self.data_hmmfile.split("\n"):
			if line.startswith("#") or line is "":
				continue
			line = split("\s+", line, 22)  # just 22 because the last can contain \s+ characters
			matrix.append(line)
			#print line[22]
		#print matrix[1][22]
		return matrix


	# Module test
if __name__ == "__main__":
	usage = """%(prog)s reads .domtblout file and returns a custom filtred \
		result (by Evalue, Bitscore, Coverage of HMM model) or it can \
		show only the Best Bitscore domain for each query on hmmscan \
		program output"""
	
	parser = argparse.ArgumentParser(description=usage)
	parser.add_argument("-b", "--bits", dest="bits",
			help="The minimum Bitscore threshold to considerate a \
			domain as a truly hits [be careful with false \
			positives] (recomended: 50)", type=float)
	parser.add_argument("-e", "--evalue", dest="evalue",
			help="The maximun E-value threshold to considerate a \
			domain as a truly hits [be careful with false \
			positives] (recomended: 1e-18)", type=float)
	parser.add_argument("-c", "--cov", dest="cov",
			help="The minimum coverage threshold to considerate a \
			domain as a truly hits [be careful with false \
			positives] (recomended: 0.35)", type=float)
	parser.add_argument("-u", "--unique", dest="unique", action="store_true",
			help="show only the Best Bitscore domain for each \
			query on hmmscan program output", default=False)
	parser.add_argument("-s", "--spacer", dest="spacer", help="Select the \
		spacer for the program output", default="\t")
	parser.add_argument("-o", "--outfile", nargs="?", help="(default: stdout)",type=argparse.FileType("w"), default=sys.stdout)

	parser.add_argument("domtblout", type=str)

	args = parser.parse_args()
	bits = args.bits
	evalue = args.evalue
	cov = args.cov
	outfile = args.outfile
	domtblout = args.domtblout
	spacer = args.spacer
	hmm = hmmParser(domtblout)

	if args.bits:
		hmm.filterByBitscore(bits)
	if args.evalue:
		hmm.filterByEvalue(evalue)
	if args.cov:
		hmm.filterByCoverage(cov)
	if args.unique:
		hmm.uniqueByBestBitscore()
	
	header = ["target name", "accession", "tlen", "query name",
		  "accession", "qlen", "E-value", "score", "bias", "#", "of",
       	          "c-Evalue", "i-Evalue", "score", "bias", "from", "to",
       	          "from", "to", "from", "to", "acc", "description of target"]

	outfile.write(spacer.join(header))
	outfile.write('\n')
	for row in hmm.matrix:
		row = spacer.join(row)
		outfile.write(row)
		outfile.write('\n')
	outfile.close()
