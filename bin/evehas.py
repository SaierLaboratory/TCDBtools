#!/usr/bin/env python

import os
import sys
import json
from Bio import AlignIO, SeqIO
import argparse
import numpy as np

import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.cm as cm
HYDROPATHY = {'G':-0.400, 
	'I':4.500, 
	'S':-0.800, 
	'Q':-3.500, 
	'E':-3.500, 
	'A':1.800, 
	'M':1.900, 
	'T':-0.700, 
	'Y':-1.300, 
	'H':-3.200, 
	'V':4.200, 
	'F':2.800, 
	'C':2.500, 
	'W':-0.900, 
	'K':-3.900, 
	'L':3.800, 
	'P':-1.600, 
	'N':-3.500, 
	'D':-3.500, 
	'R':-4.500, 
	'U':0, 
	'B':-3.500, 
	'J':-3.500, 
}

def collapse(arr):
	newarr = arr[~np.isnan(arr)]
	nans = np.nonzero(np.isnan(arr))[0]
	return newarr, nans

def uncollapse(arr, nans):
	newarr = np.zeros(len(arr) + len(nans))

	j = 0
	for i in range(len(newarr)):
		#print(i, len(newarr), nans, len(newarr) in nans)
		#if len(newarr) in nans: newarr[i] = np.nan
		if i in nans: newarr[i] = np.nan
		else: 
			newarr[i] = arr[j]
			j += 1

	return newarr

def evehas(msa, window=19):
	kernel = np.ones(window)/window
	eveh = np.zeros((len(msa), msa.get_alignment_length()))
	mlen = len(msa)
	for seqi, seq in enumerate(msa):
		rawhydro = np.array([HYDROPATHY.get(resn, np.nan) for resn in seq.seq])

		hydro, nans = collapse(rawhydro)
		if len(hydro) < window: eveh[seqi] = np.nan * np.ones(eveh.shape[1], dtype=np.float64)
		else:
			ahydro = np.convolve(hydro, kernel, 'same')

			reahydro = uncollapse(ahydro, nans)
			eveh[seqi] = reahydro

	return eveh

lcmap = cm.get_cmap('tab10')
def plot_hydropathy(hydro, ax):
	obj = ax.plot(hydro, color=lcmap(0), label='Hydropathy')
	ax.set_xlim([0, len(hydro)])
	ax.set_ylim([-3, 3])
	return obj

def plot_histogram(residues, ax, normalize=False):
	obj = None
	total = sum(residues.values())

	resn = []
	resi = []
	count = []
	for i, k in enumerate(sorted(residues)):
		#if k == '-': continue
		resi.append(i)
		resn.append(k)
		if normalize: count.append(residues[k] / total)
		else: count.append(residues[k])

	ax.set_xlim(None)
	ax.set_xlim([0, len(residues)])
	ax.set_xticks(np.arange(len(residues))+0.5)
	ax.set_xticklabels(resn)
	obj = ax.hist(resi, weights=count, bins=np.arange(-1, len(resn)+1), color=lcmap(0))
	if normalize: ax.set_ylim([0, 1])
	else: ax.set_ylim([0, total])
		
	return obj[-1]

def shannon(obs):
	classes = len(obs)
	if '-' in obs: classes -= 1

	if not classes: return 0

	total = sum(obs.values())

	entropy = 0
	for resn in obs:
		if resn != '-': continue

		entropy += -obs[resn]/total * np.log2(obs[resn]/total)

	return entropy

def plot_entropy(msa, ax, offset=-3):

	counts = []

	for row, record in enumerate(msa):

		for col, resn in enumerate(record.seq):
			if row == 0: counts.append({})

			if resn in counts[col]: counts[col][resn] += 1
			else: counts[col][resn] = 1

	values = np.array([shannon(pos) for pos in counts])

	return ax.plot(1+np.arange(len(values)), values + offset, 'tab:green', label='Diversity')

def plot_occupancy(msa, ax, offset=-3):
	counts = []

	for row, record in enumerate(msa):
		for col, resn in enumerate(record.seq):
			if row == 0: counts.append(0)
			if resn != '-': counts[col] += 1

	values = np.array(counts) / len(msa)

	return ax.plot(1+np.arange(len(values)), values + offset, 'tab:orange', label='Occupancy')
	

def get_rescounts(msa):
	columns = {}
	for seq in msa:
		for i, resn in enumerate(seq.seq):
			if i not in columns: columns[i] = {}
			if resn not in columns[i]: columns[i][resn] = 0
			columns[i][resn] += 1
	return columns

def get_occupancies(msa):
	columns = {}
	for seq in msa:
		for i, resn in enumerate(seq.seq):
			if i not in columns: columns[i] = {'gap':0, 'res':0}

			if resn in '-': columns[i]['gap'] += 1
			else: columns[i]['res'] += 1

	occupancies = {}
	for i in columns: occupancies[i] = columns[i]['res'] / (columns[i]['gap'] + columns[i]['res'])
	return occupancies

def mark_residues(arr, ax):
	row, col = np.nonzero(arr)
	outarr = np.zeros(arr.shape + (4,))
	outarr[row,col] += np.array([0.0, 1.0, 0.0, 0.5])
	
	im = ax.imshow(outarr)
	return [im]

def paint_cross(arr, row, col, ax):
	fc = (0.5, 0.5, 0.5, 0.25)
	return [ax.axvline(x=col, color=fc),
		ax.axhline(y=row, color=fc),
		ax.add_patch(patches.Circle((col, row), radius=1, fc='#ffffff00', ec='c')),
	]

	#outarr = np.zeros(arr.shape + (4,))

	#outarr[row,:] = np.array([0.5, 0.5, 0.5, 0.25])
	#outarr[:,col] = np.array([0.5, 0.5, 0.5, 0.25])
	#outarr[row,col] = np.array([0.0, 0.0, 0.0, 0.0])
	#im = ax.imshow(outarr)
	#return [im]

def paint_mark(col, ax):
	vline = ax.axvline(x=col, color='gray', lw=0.5)
	return [vline]

class EveHAS(object):
	e2 = None
	h2 = None
	h3 = None
	mode = None
	selboxes = None
	title = None
	selcross = None
	selmark = None

	fig = None
	ax = None
	ax2 = None
	ax3 = None

	names = None
	hydro = None
	msa = None
	columns = None
	histcols = None
	occupancies = None
	seqblock = None

	alnwidth = None
	alnheight = None

	normalize = False

	aspect = 1.0

	row = None
	col = None

	context = 3

	fresh = True


	def update_pos(self, row, col):
		if self.row == row and self.col == col: return
		
		#print('{} at pos. {}: {}'.format(self.names[row], col, self.msa[row].seq[col]))

		if self.row != row: 
			self.ax2.set_title('Row: {} ({}aa)'.format(self.names[row], len([resn for resn in self.msa[row] if resn != '-'])))
			if self.h2 is not None: [x.remove() for x in self.h2]
			self.h2 = plot_hydropathy(self.hydro[row], self.ax2)

			if self.fresh:
				self.ax2.legend()
				self.fresh = False

		if self.col != col:
			self.ax3.set_title('Column: {} ({:0.0%}occ)'.format(col, self.occupancies[col]))
			if self.h3 is not None: [x.remove() for x in self.h3]
			self.h3 = plot_histogram(self.columns[col], self.ax3, normalize=self.normalize)
			self.histcols = sorted(self.columns[col])

		if self.selcross is not None: [x.remove() for x in self.selcross]
		self.selcross = paint_cross(self.hydro, row, col, self.ax)

		if self.selmark is not None: [x.remove() for x in self.selmark]
		self.selmark = paint_mark(col, self.ax2)

		self.row = row
		self.col = col

		self.fig.canvas.draw()

	def do_stuff(self, event):
		if event.xdata is None or event.ydata is None: return

		if (0 <= event.xdata <= self.alnwidth) and (0 <= event.ydata <= self.alnheight):
			if event.inaxes == self.ax: 
				row = int(event.ydata + 0.5)
				col = int(event.xdata + 0.5)
				self.update_pos(row, col)

			elif event.inaxes == self.ax2:
				col = int(event.xdata + 0.5)
				if col != self.col: 
					row = self.row
					self.update_pos(row, col)

			elif event.inaxes == self.ax3:
				row = self.row
				col = int(event.xdata)

				if self.histcols is None: return
				if not (0 <= col <= len(self.histcols)-1): return

				letters = self.histcols[:]
				if self.selboxes is not None: [x.remove() for x in self.selboxes]
				self.selboxes = mark_residues(self.seqblock == letters[col], self.ax)
				self.ax.set_title('Marked {} {} residues'.format(np.sum(self.seqblock == letters[col]), letters[col]))
				self.fig.canvas.draw()
			else: return

	def show_interactive(self, msa, hydro, title=None, basename=None):
		self.msa = msa
		self.names = [seq.description for seq in self.msa]
		self.hydro = hydro
		self.title = title
		self.basename = basename

		self.alnwidth, self.alnheight = self.msa.get_alignment_length(), len(self.msa)
		#columns = np.array([[seq.seq[i] for seq in msa] for i in range(alnwidth)])

		self.columns = get_rescounts(self.msa)
		self.occupancies = get_occupancies(self.msa)
		self.seqblock = np.vstack([list(seq.seq) for seq in self.msa])

		for k in plt.rcParams:
			if k.startswith('keymap'):
				removeme = []
				for key in plt.rcParams[k]:
					if key == key.upper(): removeme.append(key)
				for key in removeme: plt.rcParams[k].remove(key)


		self.fig = plt.figure()
		self.fig.set_tight_layout(True)
		self.fig.canvas.set_window_title('EveH ({}x{}): {}'.format(self.alnheight, self.alnwidth, self.basename))

		self.ax = self.fig.add_subplot(121, aspect=self.aspect)
		self.ax2 = self.fig.add_subplot(222)
		self.ax3 = self.fig.add_subplot(224)
		self.ax3.set_ylim([0, len(self.msa)])

		self.ax.set_title(self.title)

		def onmousemove(event):
			if self.mode != 'msadrag': return

			self.do_stuff(event)

		def onmousedown(event):
			if event.inaxes == self.ax and event.button == 1: self.mode = 'msadrag'
			elif event.inaxes == self.ax2 and event.button == 1: self.mode = 'msadrag'

		def onmouseup(event):
			if self.mode == 'msadrag': self.mode = None

			self.do_stuff(event)


		def onkeyup(event):
			if event.key == None: return
			elif event.key == '?': self.print_help()
			elif event.key in 'ACDEFGHIKLMNPQRSTVWY-': 
				if self.selboxes is not None: [x.remove() for x in self.selboxes]
				self.selboxes = mark_residues(self.seqblock == event.key, self.ax)
				self.ax.set_title('Marked {} {} residues'.format(np.sum(self.seqblock == event.key), event.key))
				self.ax.set_aspect(self.aspect)
				self.fig.canvas.draw()
			elif event.key == 'c':
				if event.inaxes is self.ax:
					row = int(event.ydata)
					col = int(event.xdata)
					rowslice = slice(max(row - self.context, 0), row + self.context + 1)
					colslice = slice(max(col - self.context, 0), col + self.context + 1)
					out = 'Context:'
					for seq in self.msa[rowslice]:
						out += '\n{}'.format(seq[colslice].seq)
				elif event.inaxes is self.ax2:
					if not self.row: return
					col = int(event.xdata)
					colslice = slice(max(col - self.context, 0), col + self.context + 1)
					out = '{}'.format(self.msa[self.row][colslice].seq)
				elif event.inaxes is self.ax3:
					print('This feature is not implemented for this subplot', file=sys.stderr)
				else: return
				print(out)
			elif event.key == 'd':
				if event.inaxes is self.ax2:
					if not self.row: return
					SeqIO.write([self.msa[self.row]], sys.stdout, "fasta")
				elif event.inaxes is self.ax3:
					if not self.col: return
					print(json.dumps(self.columns[self.col]))
				elif event.inaxes is self.ax:
					print('This feature is not implemented for this subplot', file=sys.stderr)
			elif event.key == 'n': 
				self.normalize = not self.normalize
				if self.normalize: self.ax3.set_ylim([0, 1.0])
				else: self.ax3.set_ylim([0, len(self.msa)])

				if self.col is not None:
					self.ax3.set_title('Column: {} ({:0.0%}occ)'.format(self.col, self.occupancies[self.col]))
					if self.h3 is not None: [x.remove() for x in self.h3]
					self.h3 = plot_histogram(self.columns[self.col], self.ax3, normalize=self.normalize)
					self.histcols = sorted(self.columns[self.col])
				self.fig.canvas.draw()
			elif event.key == 'x':
				if event.inaxes is self.ax: 
					out = '{}\t{}'.format(self.msa[int(event.ydata)].id, self.msa[int(event.ydata)][int(event.xdata)])
				elif event.inaxes is self.ax2:
					if self.row is None: return
					out = '{}\t{}'.format(self.msa[self.row].id, self.msa[self.row][int(event.xdata)])
				elif event.inaxes is self.ax3:
					if self.col is None: return
					colrescount, colresfreq, colresnongapfreq, colresname = self.get_colrescount(int(event.xdata))
					out = '{}\t{}\t{:0.1%}\t{:0.1%}'.format(colresname, colrescount, colresfreq, colresnongapfreq)
				else: return
				print(out)
			elif event.key in ' ':
				if self.selboxes is not None: 
					[x.remove() for x in self.selboxes]
					self.selboxes = None
				self.ax.set_title(self.title)
				self.fig.canvas.draw()
			elif event.key == ']': 
				self.aspect = np.clip(self.aspect * 0.9, 0.025, 40)
				self.ax.set_aspect(self.aspect)
				self.fig.canvas.draw()
			elif event.key == '[': 
				self.aspect = np.clip(self.aspect * 1.1, 0.025, 40)
				self.ax.set_aspect(self.aspect)
				self.fig.canvas.draw()
			elif event.key == '\\':
				oldaxbox = self.ax.get_position()
				totdim = lambda box: box.x1 - box.x0 + box.y1 - box.y0

				self.aspect = np.clip(self.aspect * 1.001, 0.025, 40)
				self.ax.set_aspect(self.aspect)
				preraxbox = self.ax.get_position()

				if False: pass #close enough	
				elif totdim(preraxbox) < totdim(oldaxbox):
					self.aspect = np.clip(self.aspect / 1.1, 0.025, 40)
					self.ax.set_aspect(self.aspect)
				elif totdim(preraxbox) > totdim(oldaxbox):
					self.aspect = np.clip(self.aspect * 1.1, 0.025, 40)
					self.ax.set_aspect(self.aspect)

				self.fig.canvas.draw()


				#help(self.ax)
				#print(self.ax, self.ax.get_width() * self.ax.get_height())
				#print(self.fig)
				
		self.fig.canvas.mpl_connect('button_release_event', onmouseup)
		self.fig.canvas.mpl_connect('button_press_event', onmousedown)
		self.fig.canvas.mpl_connect('motion_notify_event', onmousemove)
		self.fig.canvas.mpl_connect('key_release_event', onkeyup)
		
		divnorm = colors.DivergingNorm(vmin=-3, vcenter=0, vmax=3)
		im = self.ax.imshow(self.hydro, cmap='RdYlBu_r', norm=divnorm)
		self.fig.colorbar(im, ax=self.ax)

		self.ax2.axhline(y=0, color='k', linewidth=0.5)
		plot_entropy(self.msa, ax=self.ax2)
		plot_occupancy(self.msa, ax=self.ax2)
		self.ax2.set_ylim([-3, 3])

		print('Press [?] for EveHAS help', file=sys.stderr)

		plt.show()

	def get_colrescount(self, resi):
		if self.col is None: return 0, 0, 'X'
		column = self.columns[self.col]
		total = 0
		nongap = 0
		value = 0
		resname = 'X'
		for i, resn in enumerate(sorted(column)):
			total += column[resn]
			nongap += column[resn] if resn != '-' else 0
			if resi == i: 
				value = column[resn]
				resname = resn
		return value, value/total, value/nongap, resname

	def print_help(self):
		print('''EveHAS help:
? - show this help printout
[ - make taller one step
] - make wider one step
\\ - autoresize one step
A-Z (capital) - highlight a specific residue
(space) - clear residue highlighting
(click) - show stats for sequence/column

c - print context, i.e. sequences from -3 to +3, residues from -3 to +3, or exact column residue counts depending on the subplot selected
d - dump a sequence as a FASTA or dump all residue counts for a specific column depending on the subplot selected
n - normalize residue histograms
x - print sequence name and column number under mouse
''')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('--window', type=int, default=19)
	parser.add_argument('--title')
	args = parser.parse_args()

	msalist = AlignIO.parse(args.infile, 'clustal')
	for msa in msalist:
		hydro = evehas(msa, window=args.window)
		evehas = EveHAS()
		evehas.show_interactive(msa, hydro, title=args.title, basename=os.path.basename(args.infile))
