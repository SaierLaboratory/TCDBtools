#!/usr/bin/env python

import os
import re
import argparse

def simplify_bn(bn, delim='_vs_', extension='.hhr'):
	newbn = bn
	newbn = newbn.replace(extension, '')
	newbn = tuple(newbn.split(delim))

	return newbn

def autoparse(val):
	try: return int(val)
	except ValueError:
		try: return float(val)
		except ValueError: return val

def purge_empties(table):
	allfams = set()
	for fam1, fam2 in table:
		allfams.add(fam1)
		allfams.add(fam2)

	onlyempty = {}
	for fam in allfams: onlyempty[fam] = True
	for fam1, fam2 in table:
		if table[fam1,fam2]:
			onlyempty[fam1] = False
			onlyempty[fam2] = False

	for fam in onlyempty:
		if onlyempty[fam]: 
			for fam1, fam2 in table:
				if fam == fam1: table.pop((fam1,fam2))
				elif fam == fam2: table.pop((fam1,fam2))
	return table

def main(indir, outdir, delim='_vs_', purge_empty=False, nomatrix=False, extension='.hhr'):
	if not os.path.isdir(outdir): os.mkdir(outdir)
	table = {}
	keys = set()
	for bn in os.listdir(indir):
		if bn.endswith(extension): 
			row = simplify_bn(bn, delim)
			table[row] = {}
			with open('{}/{}'.format(indir, bn)) as fh:
				for l in fh:
					if l.startswith('Match_columns'):
						k, v = l.strip().split()
						k = 'Query_cols'
						if k not in table[row]: table[row][k] = autoparse(v)
						keys.add(k)

					elif l.startswith('  1'):
						k = 'Template_cols'
						v = l.strip().split()[-1][1:-1]
						if k not in table[row]: table[row][k] = autoparse(v)
						keys.add(k)
					elif l.startswith('Probab='):
						for kv in l.strip().split():
							k, v = kv.split('=')
							if k not in table[row]: table[row][k] = autoparse(v)
							keys.add(k)
	for row in table:
		table[row]['Query_cov'] = '{:0.1%}'.format(table[row]['Aligned_cols'] / table[row]['Query_cols'])
		table[row]['Template_cov'] = '{:0.1%}'.format(table[row]['Aligned_cols'] / table[row]['Template_cols'])
		keys.add('Query_cov')
		keys.add('Template_cov')
			
	if purge_empty: purge_empties(table)

	allfams = set()
	for fam1, fam2 in table: 
		allfams.add(fam1)
		allfams.add(fam2)

	for k in keys:
		if not nomatrix:
			with open('{}/{}_matrix.tsv'.format(outdir, k), 'w') as fh:
				for fam1 in sorted(allfams): fh.write('\t{}'.format(fam1))
				fh.write('\n')

				for fam1 in sorted(allfams):
					fh.write(str(fam1))
					for fam2 in sorted(allfams):
						cell = table.get((fam1,fam2), {})
						value = cell.get(k, '')
						fh.write('\t{}'.format(value))
					fh.write('\n')
		with open('{}/{}_columnwise.tsv'.format(outdir, k), 'w') as fh:
			for fam1 in sorted(allfams):
				for fam2 in sorted(allfams):
					if (fam1,fam2) in table:
						fh.write('{}\t{}\t{}\n'.format(fam1, fam2, table[fam1,fam2].get(k, '')))


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('indir', default='.', nargs='?', help='Directory with .hhr results')
	parser.add_argument('-d', '--delimiter', default='_vs_', help='Delimiter')
	parser.add_argument('-o', '--outdir', required=True, help='Output directory')
	parser.add_argument('-u', action='store_true', help='Remove empty rows/columns')
	parser.add_argument('--no-matrix', action='store_true', help='Skip reshaping columnwise table into a matrix')
	parser.add_argument('-x', '--extension', default='.hhr', help='Extension for HHalign results')

	args = parser.parse_args()

	main(args.indir, args.outdir, delim=args.delimiter, purge_empty=args.u, nomatrix=args.no_matrix, extension=args.extension)
