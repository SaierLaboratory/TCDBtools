#!/usr/bin/env python

from __future__ import print_function, division
import os
if 'MPLBACKEND' not in os.environ: os.environ['MPLBACKEND'] = 'Qt4Agg'

import argparse
import quod
import numpy as np

import sys
import subprocess
import re
import shlex

#TODO: generalize key classes and put them together
import phoboshop

#TODO: somehow make this forward-compatible
import urllib

import warnings
try: import avehas3
except ImportError: 
	avehas3 = None
	warnings.warn('Could not find avehas3. MSA plotting disabled')

def TMCOLOR(i):
	if i is None: return 'orange'

	elif i % 3 == 0: return 'orange'
	elif i % 3 == 1: return 'cyan'
	elif i % 3 == 2: return 'purple'

	else: return 'orange'

def LINECOLOR(i):
	if i is None: return 'red'

	elif i % 3 == 0: return 'red'
	elif i % 3 == 1: return 'blue'
	elif i % 3 == 2: return 'green'

	else: return 'red'

class Protein(phoboshop.Protein):
	pass
	vspans = []
	what = None
	spans = []

	def get_spans(self):
		#hmmtop = self.what.entities[1]
		#return hmmtop.spans
		return self.spans

	def plot_vspans(self, ax, style=0):
		color = TMCOLOR(style)
		alpha = 0.25

		vspans = []
		for span in self.spans:
			vspans.append(ax.axvspan(span[0], span[1], facecolor=color, alpha=alpha))
		self.vspans = vspans

	def render(self, style=0):
		self.what = quod.What(self.get_sequence(), style=style)
		self.entities['what'] = self.what
		entlist = []
		entlist = [self.what.entities[0]]

		return entlist

	def hmmtop(self):
		fasta = '>untitled\n{}'.format(self.get_sequence())
		p = subprocess.Popen(['hmmtop', '-if=--', '-pi=spred', '-sf=FAS', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		out, err = p.communicate(input=fasta)


		corrections = []
		cumul = 0
		for resn in self.get_sequence():
			if resn not in 'ACDEFGHIKLMNPQRSTVWY': cumul += 1
			corrections.append(cumul)

		if not out.strip(): self.spans = []
		else:
			s = re.findall('(?:(?:[0-9]+ *)+)$', out.strip())[0]
			rawindices = [int(x) for x in s.split()[1:]]
			indices = [[rawindices[i], rawindices[i+1]] for i in range(0, len(rawindices), 2)]
			self.spans = indices

		for tms in self.spans:
			for i, index in enumerate(tms):
				tms[i] += corrections[index-1]

	def get_ydro(self, targetx):
		hydro = self.what.entities[0]
		for x, y in zip(hydro.X, hydro.Y):
			if x == targetx: return y
		return 0

class Alignment(Protein):
	vspans = []
	what = None
	spans = []
	alignments = None
	header = 'alignment'
	hydro = None
	amphi = None
	tmcenters = None
	simil = None
	amphipathicity = False
	def __init__(self, alignments):
		self.alignments= alignments

	def get_alignments(self): return self.alignments

	def render(self, style=0):
		for msa in self.alignments:
			self.what = quod.What('', style=style, nohmmtop=True)
			self.hydro = self.what.entities[0]
			self.hydro.Y = avehas3.get_average_hydropathies(msa, window=19)
			self.hydro.X = np.arange(0., len(self.hydro.Y))

			self.entities['what'] = self.what
			entlist = []
			entlist = [self.what.entities[0]]

			if self.amphipathicity:
				amphi = quod.Hydropathy('', style=style+2)
				amphi.Y = avehas3.get_average_amphipathicities(msa, window=19)
				amphi.X = np.arange(0, len(amphi.Y))
				entlist.append(amphi)

			entlist.append(avehas3.TMcenter(avehas3.get_tmcenters(msa), ymin=-3, ymax=0.1))

			simil = avehas3.get_similarities(msa)
			window = 10
			ymin = -3
			ymax = -1.5
			scores = np.array([np.nanmean(simil[i:i+window]) for i in range(0, len(simil)-window)])
			normscores = (scores - min(scores)) / (max(scores) - min(scores))
			scaledscores = normscores * (ymax - ymin)
			shiftedscores = scaledscores + ymin
			similcurve = quod.Hydropathy('', style='gray')
			similcurve.Y = shiftedscores
			similcurve.X = np.arange(0, len(shiftedscores)) + window//2
			entlist.append(similcurve)
			

			#what's the worst that can happen???
			return entlist


class FragmentProtein(Protein):
	def __init__(self, fragfasta, fullfasta):
		if not fragfasta.startswith('>'): raise ValueError('fasta is not in FASTA format')
		elif '\n' not in fragfasta: raise ValueError('fasta is not in FASTA format')
		if not fullfasta.startswith('>'): raise ValueError('fasta is not in FASTA format')
		elif '\n' not in fullfasta: raise ValueError('fasta is not in FASTA format')

		self.fasta = fragfasta
		self.fragfasta = fragfasta
		self.fullfasta = fullfasta
		self.header = self.get_header()
		self.sequence = self.get_sequence()

	def render(self, style=0):
		self.what = quod.FragmentWhat(self.fragfasta, self.fullfasta, style=style)
		self.entities['what'] = self.what
		entlist = []
		entlist = [self.what.entities[0]]

		return entlist

	def hmmtop(self):
		fasta = self.fragfasta

		p = subprocess.Popen(['hmmtop', '-if=--', '-pi=spred', '-sf=FAS', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		out, err = p.communicate(input=fasta)

		corrections = []
		cumul = 0
		for resn in self.get_sequence():
			if resn not in 'ACDEFGHIKLMNPQRSTVWY': cumul += 1
			corrections.append(cumul)

		if not out.strip(): self.spans = []
		else:
			s = re.findall('(?:(?:[0-9]+ *)+)$', out.strip())[0]
			rawindices = [int(x) for x in s.split()[1:]]
			indices = [[rawindices[i], rawindices[i+1]] for i in range(0, len(rawindices), 2)]
			self.spans = indices

		for tms in self.spans:
			for i, index in enumerate(tms):
				tms[i] += corrections[index-1]

class TMWeaver(object):
	proteins = []
	tmss = []
	mode = 'normal'
	outfmt = 'multiquod'
	fig = None
	plot = None
	pid = None
	target = None

	outfile = '/dev/stdout'
	pipeto = None
	append = False

	testline = None
	selections = {}
	queue = {}

	def __init__(self):
		pass

	def add_protein(self, protein):
		self.proteins.append(protein)
		self.tmss.append(protein.hmmtop())

	def add_alignment(self, alignment):
		self.proteins.append(alignment)

	def run(self):
		print('Use [?] to get help on phoboshop shortcuts', file=sys.stderr)
		self.fig = quod.plt.figure()
		self.plot = quod.Plot(fig=self.fig)

		for i, p in enumerate(self.proteins):
			for e in p.render(style=i):
				self.plot.add(e)
			p.plot_vspans(self.plot.ax, style=i)

		self.update_title()

		self.plot.ax.figure.canvas.mpl_connect('button_press_event', self.onmousedown)
		self.plot.ax.figure.canvas.mpl_connect('button_release_event', self.onmouseup)
		self.plot.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmousemove)
		self.plot.ax.figure.canvas.mpl_connect('key_press_event', self.onkeydown)

		self.testline, = self.plot.ax.plot([], [])
		self.plot.render()
		quod.plt.show()

	def _extract_numeric(self, string):
		nums = re.findall('[0-9]+', string)
		if nums is None: return None
		else: return int(nums[0])

	def _is_in_span(self, x, pid=0):
		inside = []
		if pid is None:
			for i, p in enumerate(self.proteins):
				for s, span in enumerate(p.get_spans()):
					if span[0] <= x <= span[1]:
						inside.append((i, s))
		else:
			for s, span in enumerate(self.proteins[pid].get_spans()):
				if span[0] <= x <= span[1]: 
					inside.append((pid, s))
		#return inside
		try: return inside[0]
		except IndexError: return None

	def _is_on_span_edge(self, x, pid=0, tolerance=2):
		onedge = []
		if pid is None:
			for i, p in enumerate(self.proteins):
				for s, span in enumerate(p.get_spans()):
					if span[0]-tolerance <= x <= span[0]+tolerance:
						onedge.append((i, s, 0))
					if span[1]-tolerance <= x <= span[1]+tolerance:
						onedge.append((i, s, 1))
		else:
			for s, span in enumerate(self.proteins[pid].get_spans()):
				if span[0]-tolerance <= x <= span[0]+tolerance:
					onedge.append((pid, s, 0))
				if span[1]-tolerance <= x <= span[1]+tolerance:
					onedge.append((pid, s, 1))

		#return onedge
		try: return onedge[0]
		except IndexError: return None

	def onmousedown(self, event): 
		if not event.inaxes: return

		if self.mode.startswith('edit'):
			#TODO: decide what to do/add onto if target is None and there are multiple proteins
			target = 0

			onedge = self._is_on_span_edge(event.xdata, self.pid)
			if onedge:
				self.mode += 'resizedrag'
				self.target = onedge
				return

			inspan = self._is_in_span(event.xdata, self.pid)
			if inspan:
				self.mode += 'movedrag'
				self.target = tuple(list(inspan) + [event.xdata])
				return

			self.mode += 'createdrag'
			vspan = self.plot.ax.axvspan(event.xdata, event.xdata+1, fc=TMCOLOR(self.pid), alpha=0.25)
			if self.pid is None: 
				self.proteins[0].vspans.append(vspan)
				self.target = (0, len(self.proteins[0].vspans)-1)
			else: 
				self.proteins[self.pid].vspans.append(vspan)
				self.target = (self.pid, len(self.proteins[self.pid].vspans)-1)
			self.plot.ax.figure.canvas.draw()

		elif self.mode.startswith('mark'):
			self.mode += 'drag'
			target = []
			for i, p in enumerate(self.proteins):
				x = int(round(event.xdata))
				y = p.get_ydro(x)
				marker = self.plot.ax.scatter([x], [y], c='k')
				target.append(marker)
				try: resn = p.get_sequence()[x]
				except IndexError: resn = ''
				label = self.plot.ax.text(x, y, resn)
				target.append(label)
			self.target = target
			
			self.plot.ax.figure.canvas.draw()

	def onmouseup(self, event): 
		#if 'drag' not in self.mode: return
		if 'resizedrag' in self.mode:
			self.mode = self._cleave(self.mode, 'resizedrag')
			vspan = self.proteins[self.target[0]].vspans[self.target[1]]
			xy = vspan.get_xy()
			x1 = xy[1,0]
			x2 = xy[2,0]
			xmin = int(round(min(x1, x2)))
			xmax = int(round(max(x1, x2)))

			self.proteins[self.target[0]].spans[self.target[1]] = [xmin, xmax]
			
			xy[:2,0] = min(x1, x2)
			xy[2:4,0] = max(x1, x2)
			xy[4:,0] = min(x1, x2)
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()

		elif 'movedrag' in self.mode:
			self.mode = self._cleave(self.mode, 'movedrag')
			vspan = self.proteins[self.target[0]].vspans[self.target[1]]
			xy = vspan.get_xy()
			x1 = xy[1,0]
			x2 = xy[2,0]

			xmin = int(round(min(x1, x2)))
			xmax = int(round(max(x1, x2)))

			self.proteins[self.target[0]].spans[self.target[1]] = [xmin, xmax]
			vspan.figure.canvas.draw()

		elif 'createdrag' in self.mode:
			self.mode = self._cleave(self.mode, 'createdrag')
			vspan = self.proteins[self.target[0]].vspans[self.target[1]]
			xy = vspan.get_xy()
			x1 = xy[1,0]
			x2 = xy[2,0]

			xmin = int(round(min(x1, x2)))
			xmax = int(round(max(x1, x2)))
			xy[:2,0] = min(x1, x2)
			xy[2:4,0] = max(x1, x2)
			xy[4:,0] = min(x1, x2)
			vspan.set_xy(xy)

			self.proteins[self.target[0]].spans.append([xmin, xmax])
			vspan.figure.canvas.draw()

		elif 'markdrag' in self.mode:
			self.mode = self._cleave(self.mode, 'drag')
			n = 0
			for i in range(0, len(self.target), 2):
				marker = self.target[i]
				label = self.target[i+1]
				print('#{}{:d} in {} (= {:0.2f}kcal/mol)'.format(label.get_text(), int(round(marker.get_offsets()[0,0])), self.proteins[n].header, marker.get_offsets()[0,1]), file=sys.stderr)
				marker.remove()
				label.remove()
				n += 1
			self.target = None
			self.plot.ax.figure.canvas.draw()

		elif 'delete' in self.mode:
			inspan = self._is_in_span(event.xdata, pid=None)
			if not inspan: return

			vspan = self.proteins[inspan[0]].vspans.pop(inspan[1])
			self.proteins[inspan[0]].spans.pop(inspan[1])
			vspan.remove()
			self.plot.ax.figure.canvas.draw()
			
	def _cleave(self, string, suffix):
		if suffix in string:
			return string[:string.find(suffix)]
		else: return string

	def onmousemove(self, event): 
		if 'drag' not in self.mode: return
		elif event.inaxes != self.plot.ax: return

		elif 'resizedrag' in self.mode:
			vspan = self.proteins[self.target[0]].vspans[self.target[1]]
			side = self.target[2]
			xy = vspan.get_xy()
			if side == 0:
				xy[:2,0] = event.xdata
				xy[4:,0] = event.xdata
			elif side == 1:
				xy[2:4,0] = event.xdata
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()

		elif 'movedrag' in self.mode:
			vspan = self.proteins[self.target[0]].vspans[self.target[1]]
			offset = event.xdata - self.target[2]
			xy = vspan.get_xy()
			xy[:,0] += offset
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()

			self.target = tuple(list(self.target[:2]) + [event.xdata])
			
		elif 'createdrag' in self.mode:
			vspan = self.proteins[self.target[0]].vspans[self.target[1]]
			xy = vspan.get_xy()
			xy[2:4,0] = event.xdata
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()

		elif 'markdrag' in self.mode:
			n = 0
			#for i, marker in enumerate(self.target):
			#for i, marker in enumerate(self.target):
			for i in range(0, len(self.target), 2):
				x = int(round(event.xdata))
				y = self.proteins[n].get_ydro(x)
				marker = self.target[i]
				label = self.target[i+1]
				marker.set_offsets([[x,y]])

				try: resn = self.proteins[n].get_sequence()[x]
				except IndexError: resn = ''
				label.set_position([x, y])
				label.set_text(resn)

				n += 1
			self.plot.ax.figure.canvas.draw()

	def onkeydown(self, event): 
		if 'drag' in self.mode: return

		if event.key == '?':
			self.print_help()

		elif re.match('[0-9]+', event.key):
			#TODO: generalize for >10 sequences
			pid = int(event.key)
			if pid > len(self.proteins): pid = len(self.proteins)
			self.mode = 'edit{}'.format(pid)
			self.pid = pid - 1
		elif event.key == 'a':
			self.mode = 'edit'
			self.pid = None
		elif event.key == 'i':
			self.mode = 'edit'
			self.pid = None

		elif event.key == 'd':
			self.mode = 'delete'

		elif event.key == 'm':
			self.mode = 'mark'

		elif event.key == 'w':
			self.write_tmss()

		elif event.key == 'escape':
			self.mode = 'normal'

		self.update_title()

	def write_tmss(self):
		out = ''
		allseq = ''
		for protein in self.proteins:
			allseq += protein.fasta + '\n'
			if self.outfmt == 'fasta':
				s = ''
				for i, span in enumerate(sorted(protein.get_spans())):
					s += '{}_TMS{}\n'.format(protein.header, i+1)
					s += protein.get_sequence()[span[0]-1:span[1]]
					s += '\n'
				out += s + '\n'
			elif self.outfmt == 'hmmtop':
				s = '>HP: {} {}  UNK  {}  '.format(len(protein.get_sequence()), protein.header, len(protein.get_spans()))
				for span in sorted(protein.get_spans()):
					s += '{}  {}  '.format(span[0], span[1])
				out += '{}\n'.format(s)
			elif self.outfmt == 'multiquod':
				s = '#add TMSs for {}\n'.format(protein.header)
				s += 'tms add SUBPLOT COLOR '
				for span in sorted(protein.get_spans()):
					s += ' {} {}'.format(span[0], span[1])
				out += '{}\n'.format(s)

			###!!!
			elif self.outfmt == 'ranges':
				s = ''
				for span in sorted(protein.get_spans()):
					s += '{}-{},'.format(span[0], span[1])
				if s.endswith(','): s = s[:-1]
				out += s + '\n'
			elif self.outfmt == 'tsv':
				s = ''
				for span in sorted(protein.get_spans()):
					for x in span: s += '{}\t'.format(x)
				out += s.strip() + '\n'

		obj = {'user': self.user,
			'seq': protein.get_sequence(),
			'tmss': out.strip()
		}
		f = urllib.urlopen(self.outfile, data=urllib.urlencode(obj))
		resp = f.read()
		print(resp)

	def update_title(self):
		title = 'Insert title here'
		if self.mode == 'normal': title = 'NORMAL mode'
		elif self.mode.startswith('edit'): 
			if self.pid is None: title = 'EDIT mode'
			else: title = 'EDIT {} mode'.format(self.pid+1)
		elif self.mode.startswith('delete'): title = 'DELETE mode'
		elif self.mode.startswith('mark'): title = 'MARK mode'

		self.fig.canvas.set_window_title(title)

	def print_help(self):
		print('''
CONTROLS
========

1-9 - Enter edit mode for a specific sequence
a - Enter edit mode
d - Enter deletion mode
i - Enter edit mode
m - Enter mark mode
w - Write cut sequence(s) to disk/stdout. TODO: run the selected pipeto program if defined (useful for TMSOC, presumably)
? - Pull up this help page
ESC - Enter normal mode, which does absolutely nothing
	''', file=sys.stderr)

def detect_filetype(f):
	f.seek(0)
	for l in f:
		if not l.strip(): continue
		elif l.startswith('CLUSTAL'): return 'clustal'
		elif l.startswith('>'): return 'fasta'
		else: return 'rawseq'

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('username')
	parser.add_argument('infile', nargs='+', help='Sequence files to read')
	if 'TMTRAINER' in os.environ: parser.add_argument('-o', metavar='ENDPOINT', default=os.environ['TMTRAINER'], help='URL pointing to a TMTrainer receptacle (default: $TMTRAINER == {})'.format(os.environ['TMTRAINER']))
	else: parser.add_argument('-o', metavar='ENDPOINT', help='URL pointing to a TMTrainer receptacle. Set $TMTRAINER to avoid having to specify this on each run.')
	parser.add_argument('--amphipathicity', action='store_true', help='Also plot amphipathicity (currently only implemented for MSAs)')
	#parser.add_argument('--frag', action='store_true', help='Open sequences as fragments (frag1 full1 frag2 full2... fragN fullN)')
	#TODO: decide how to handle overlapping "TMSs"
	args = parser.parse_args()

	if not args.o:
		print('Please specify an endpoint')
		exit(1)

	tmweaver = TMWeaver()
	tmweaver.outfile = args.o
	#tmweaver.append = args.a

	tmweaver.outfmt = 'ranges'
	tmweaver.user = args.username

	#if args.frag:
	if False:
		for i in range(0, len(args.infile), 2):
			fragfn = args.infile[i]
			if fragfn.startswith('asis:'): fragseq = fragfn[5:]
			else: 
				with open(fragfn) as f: fragseq = f.read()
			if not fragseq.startswith('>'): fragseq = '>fragseq{}\n{}'.format(i+1, fragseq)

			fullfn = args.infile[i+1]
			if fullfn.startswith('asis:'): fullseq = fullfn[5:]
			else: 
				with open(fullfn) as f: fullseq = f.read()
			if not fullseq.startswith('>'): fullseq = '>fullseq{}\n{}'.format(i+1, fullseq)

			p = FragmentProtein(fragseq, fullseq)
			tmweaver.add_protein(p)
	else:
		for infn in args.infile:
			f = open(infn)
			filetype = detect_filetype(f)
			f.seek(0)
			if filetype == 'fasta':
				for fasta in phoboshop.split_fasta(f):
					p = Protein(fasta=fasta)
					tmweaver.add_protein(p)
				f.close()
			elif filetype == 'clustal':
				if avehas3 is None: raise ImportError('Could not find AveHAS3')
				alignments = avehas3.AlignIO.parse(f, 'clustal')
				aln = Alignment(alignments=alignments)
				aln.amphipathicity = args.amphipathicity
				aln.header = infn
				tmweaver.add_alignment(aln)
					
					
	tmweaver.run()
