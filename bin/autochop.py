#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse
import scipy.signal
import numpy as np
import Bio.SeqIO as SeqIO
import Bio.Seq as Seq
import warnings
warnings.filterwarnings('ignore')
from quod import HYDROPATHY
from kevtools_common import badnan

def get_long_nonzeros(trutharr, length=1):
	suffic = np.ones(trutharr.shape)
	for i in range(int(length)):
		suffic *= np.roll(trutharr, -i)
	results = np.zeros(trutharr.shape)
	for i in range(int(length)):
		results += np.roll(suffic, i)
	results = np.clip(results, 0, 1)

	return results

def autochop(seq, **kwargs):

	window = kwargs.get('window', 15)
	min_length = kwargs.get('min_length', 30)
	max_hydro = kwargs.get('max_hydro', -1)
	interactive = kwargs.get('interactive', False)

	rawhydro = np.array([HYDROPATHY.get(resn.upper(), np.nan) for resn in seq.seq])
	nrawhydro, nans = badnan.collapse(rawhydro)
	b, a = scipy.signal.butter(3, 1/window, btype='lowpass', output='ba')
	nfilthydro = scipy.signal.lfilter(b, a, nrawhydro)
	filthydro = badnan.uncollapse(nfilthydro, nans)

	deleteme = get_long_nonzeros(filthydro < max_hydro, length=args.min_length)

	newaa = ''.join([seq.seq[i] for i in range(len(seq.seq)) if not deleteme[i]])

	newseq = SeqIO.SeqRecord(seq=Seq.Seq(newaa), id=seq.id, name=seq.name, description=seq.description)

	if interactive:
		kernel = np.ones(int(round(window)))/int(round(window))
		
		navehydro = np.convolve(nrawhydro, kernel, 'same')
		avehydro = np.roll(badnan.uncollapse(navehydro, nans), int(round(window)/2)+4)
		plt.plot(filthydro, color='gray', alpha=0.5)
		plt.plot(avehydro)
		plt.plot(avehydro / deleteme, color='red')
		plt.ylim([-3,3])
		plt.axhline(y=0, lw=1, color='k')
		plt.show()
	
	return newseq

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='*', default=['/dev/stdin'], help='FASTAs to read in (default:stdin)')
	parser.add_argument('--window', type=float, default=15, help='-6dB cutoff for low pass (default:15)')
	parser.add_argument('--min-length', type=float, default=30, help='Minimum length of contiguous segments below hydropathy cutoff')
	parser.add_argument('--max-hydro', type=float, default=-1, help='Minimum hydropathy permitted (=maximum hydropathy to chop away) (default:-1)')
	parser.add_argument('-o', default='/dev/stdout', help='Where to write the resulting sequence(s) (default:stdout)')

	parser.add_argument('--show', action='store_true', help='Show off the cuts')

	args = parser.parse_args()

	outfh = open(args.o, 'w')
	for fn in args.infile:
		for seq in SeqIO.parse(fn, 'fasta'):
			newseq = autochop(seq, window=args.window, min_length=args.min_length, max_hydro=args.max_hydro, interactive=args.show)
			outfh.write(newseq.format('fasta'))

