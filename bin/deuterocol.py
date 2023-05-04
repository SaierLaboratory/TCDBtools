#!/usr/bin/env python

import libdeuterocol
import argparse
import os

def run_deuterocol(famlist1, famlist2, bundle, tmdatadir, outdir, level=3, allow_internal=False, cut_loops=0, num_threads=1, aligner='tmalign'):
	family1 = libdeuterocol.deuterocol1.load_family(famlist1, tmdatadir, bundle=bundle, level=level)
	#family1.write('{}/family1'.format(outdir))
	family2 = libdeuterocol.deuterocol1.load_family(famlist2, tmdatadir, bundle=bundle, level=level)
	#family2.write('{}/family2'.format(outdir))

	if not os.path.isdir(outdir): os.mkdir(outdir)

	libdeuterocol.deuterocol1.write([family1, family2], '{}/deuterocol1'.format(outdir))

	deut2 = libdeuterocol.deuterocol2.Deuterocol2([family1, family2], tmdatadir=tmdatadir, bundle=bundle, cut_loops=cut_loops)
	deut2.compute()
	deut2.write(outdir)

	trypsin = libdeuterocol.trypsin.Trypsin(deut2)
	trypsin.cut_pdbs('{}/cut_pdbs'.format(outdir))

	aligntool = libdeuterocol.aligner.get_aligner(aligner)(deut2, outdir)

	aligntool.run()
	aligntool.tabulate()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('configfile', nargs='?')
	
	#TODO: Make an argument group for reading the data from a dscript
	parser.add_argument('--fams1', nargs='+', required=True)
	parser.add_argument('--fams2', nargs='+', required=True)
	parser.add_argument('--bundle', type=int, required=True)
	parser.add_argument('--tmdatadir', required=True)
	parser.add_argument('-o', required=True)
	parser.add_argument('--level', type=int, default=3)
	parser.add_argument('--allow-internal', action='store_true')
	parser.add_argument('--cut-loops', type=int, default=0)

	args = parser.parse_args()

	run_deuterocol(args.fams1, args.fams2, args.bundle, args.tmdatadir, args.o, level=args.level, allow_internal=args.allow_internal, cut_loops=args.cut_loops)

