#!/usr/bin/env python
import argparse
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-l', help='Manual title specification')
	parser.add_argument('-q', action='store_true', help='More magic')
	parser.add_argument('-s', action='store_true', help='Deprecated. Use asis:ASDFGHIKL instead')
	parser.add_argument('--correlate', action='store_true', help='Calculate pairwise correlations')
	parser.add_argument('--grid', action='store_true', help='Enable grid')
	parser.add_argument('--kernel', help='Which smoothing kernel to use, if any (default: moving-average)')
	parser.add_argument('--xticks', type=int, default=None, help='x-tick spacing')
	parser.add_argument('--width', type=float, default=None, help='Figure width')
	parser.add_argument('--height', type=float, default=None, help='Figure height')
	parser.add_argument('-o', help='Where to save the figure')
	parser.add_argument('infiles', nargs='+', help='Input files given as FRAG1 FULL1 FRAG2 FULL2 ... . Use asis: to input sequences directly as CLI arguments')

	args = parser.parse_args()

	if args.s: raise ZeroDivisionError('-s is deprecated. Use asis:ASDFGHIKL instead')

	if len(args.infiles) % 2: 
		print('Please provide an even number of sequences')
		parser.print_usage()
		exit(1)

	if args.o: os.environ['MPLBACKEND'] = 'Agg'
	elif 'MPLBACKEND' not in os.environ: os.environ['MPLBACKEND'] = 'Qt4Agg'
	if args.q: os.environ['MPLBACKEND'] = 'Agg'
	import quod
	#import matplotlib.pyplot as plt
	whats = []
	plt = quod.plt

	sequences = []
	for i, thing in enumerate(args.infiles):
		if thing.startswith('asis:'): sequences.append('>seq{}\n{}'.format(i+1, thing[5:]))
		else:
			with open(thing) as f: 
				s = f.read()
				if not s.startswith('>'): sequences.append('>seq{}\n{}'.format(i+1, s))
				else: sequences.append(s)

	fig = plt.figure()
	plot = quod.Plot(fig=fig)
	n = 0
	for i in range(0, len(sequences), 2):
		x = quod.FragmentWhat(sequences[i], sequences[i+1], style=n, kernel=args.kernel)
		plot.add(x)
		whats.append(x)
		n += 1
	plot.render()

	if args.correlate and (len(whats) > 1):
		for i in range(len(whats)):
			hydro1 = whats[i].entdict['hydro'].Y
			for j in range(len(whats)):
				if i >= j: continue

				hydro2 = whats[j].entdict['hydro'].Y

				nogap1 = []
				nogap2 = []
				for r1, r2 in zip(hydro1, hydro2):
					if quod.np.isnan(r1) or quod.np.isnan(r2): continue
					nogap1.append(r1)
					nogap2.append(r2)

				pearson = quod.np.corrcoef(nogap1, nogap2)[1,0]

				print('Sequence {} (len. {}) vs {} (len. {}): {:0.3f}'.format(i+1, len(hydro1), j+1, len(hydro2), pearson))


	#Last minute stuff
	if args.grid: plot.ax.grid(1)

	if args.xticks:
		xlim = plot.ax.get_xlim()
		plot.ax.set_xticks(quod.np.arange(xlim[0], xlim[1], args.xticks))

	if args.width: plot.fig.set_figwidth(args.width)

	if args.height: plot.fig.set_figheight(args.height)

	if args.l: plot.ax.set_title(args.l)

	if args.o: fig.savefig(args.o)
	else: plt.show()
