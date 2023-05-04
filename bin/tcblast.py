#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import subprocess, tempfile, os, sys
import numpy as np

import re

import warnings
warnings.filterwarnings("ignore")


def blast(seq, maxhits=50):
	f = tempfile.NamedTemporaryFile(delete=False)
	f.write(seq.encode('utf-8'))
	f.close()
	try:
		cmd = ['blastp', '-db', 'tcdb', '-evalue', '1.0', '-query', f.name, '-gapopen', '11', '-gapextend', '1', '-matrix', 'BLOSUM62', '-comp_based_stats', '0', '-seg', 'no']
		tabout = subprocess.check_output(cmd + ['-outfmt', '6', '-max_target_seqs', str(maxhits)])
		pairwise = subprocess.check_output(cmd + ['-outfmt', '0', '-num_descriptions', str(maxhits), '-num_alignments', str(maxhits)])
	finally: os.remove(f.name)

	return tabout.decode('utf-8'), pairwise.decode('utf-8')

def sanitize(putative_filename, tiny=0):
	if tiny: putative_filename = putative_filename.split()[0]
	fn = re.sub('[>/]', '', putative_filename)
	fn = re.sub('[^A-Za-z0-9\.]', '_', fn)
	return fn

def hmmtop(seq, outdir=None, silent=False):
	#simplify FASTAs for less ambiguous parsing:
	oldseq = seq.split('\n')[0]
	if seq.strip().startswith('>'): seq = '>seq\n' + seq[seq.find('\n')+1:]

	if outdir:
		if not os.path.isdir(outdir): os.mkdir(outdir)
		fn = str(outdir + '/' + sanitize(oldseq.split()[0]) + '.top')
	if outdir and os.path.isfile(fn):
		f = open(fn)
		x = f.read().decode('utf-8')
		f.close()
		return x

	f = subprocess.Popen(['hmmtop', '-if=--', '-sf=FAS', '-pi=spred', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	topout, err = f.communicate(input=seq)
	print(err.strip(), file=sys.stderr)

	if outdir:
		f = open(fn, 'w')
		f.write(topout)
		f.close()

	return topout.decode('utf-8')

def parse_hmmtop(topout):
	skip = 0
	found = 0
	indices = []

	#indices = re.findall('(?:\s+([0-9]+))\s*$', topout)
	indices = re.findall('((?:[0-9]+\s+)+)$', topout)[0].strip().split()[1:]

	#for i, x in enumerate(topout.split()): 
	#	if skip:
	#		skip += 1
	#		continue
	#	if x in ('IN', 'OUT'): 
	#		try: int(topout.split()[i+1])
	#		except ValueError: continue
	#		if int(topout.split()[i+1]): 
	#			indices = topout.split()[i+2:]
	#			break
	tmss = []
	for i in range(0, len(indices), 2): 
		tmss.append([int(indices[i]), int(indices[i+1])-int(indices[i])])
	return tmss

def plot_tab(tab, top, filename, maxaln=50, dpi=100, imgfmt='png', overwrite=False, silent=False):
	bars = {}
	targets = []
	maxlen = 0
	mini = {}
	#parse BLAST output
	evalue = ''

	if not overwrite and os.path.isfile(filename): return

	for l in tab.split('\n'): 
		if len(targets) >= maxaln: break

		if l.startswith('#'): continue
		hit = l.split('\t')
		if len(hit) < 10: continue
		# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
		#addme = list(map(int, hit[6:8]))
		addme = [int(hit[6]), int(hit[3])]

		#if this is the first time this target sequence appeared, check its evalue
		if hit[1] not in bars:
			if not evalue and float(hit[10]) >= 0.05: 
				evalue = hit[1]

		#add the most recent coordinates to the list of bars to draw
		try: bars[hit[1]].append(addme)
		except KeyError: bars[hit[1]] = [addme]

		#add the target sequence to the drawing queue
		if hit[1] not in targets: targets.append(hit[1])

		#adjust the maximum length of the graph based on the final covered residue?
		#maybe work in query length instead
		if sum(addme) > maxlen: maxlen = sum(addme)

		#mini: label left anchor
		try: mini[hit[1]] = min(mini[hit[1]], int(hit[6]))
		except KeyError: mini[hit[1]] = int(hit[6])


	#len(targets): # sequences to graph

	covheight = len(targets)/12.0*2
	topheight = 0.25
	bottomheight = 1.0
	
	plt.figure(figsize=(8, covheight+topheight+bottomheight), tight_layout=1, dpi=dpi)
	#gs = gridspec.GridSpec(2, 1, height_ratios=[topheight, covheight])
	gs = gridspec.GridSpec(1, 1)
	ax1 = plt.subplot(gs[0])
	ax1.axes.get_yaxis().set_visible(0)
	ax1.xaxis.tick_top()

	yoff = 3
	for i, t in enumerate(targets): 
		ax1.broken_barh(bars[t], [-i-yoff, 1*0.8])
		plt.annotate(t, [mini[t], -i-yoff+0.1], color='white', size=8)
		if t == evalue:
			#draw the 0.05 cutoff line
			#maybe spread the word that the red line is a 0.05 cutoff line?
			ax1.axhline(y=(-i-yoff-0.15), color='red')
		#help(plt.annotate)
	plt.xlim(left=0, right=maxlen)
	plt.ylim(top=2, bottom=-i-1-yoff)
	#plt.ylim(top=2, bottom=-i-1-yoff)

	#ax1 = plt.subplot(gs[0], sharex=ax1)
	#ax1.axes.get_xaxis().set_visible(0)
	ax1.broken_barh(parse_hmmtop(top), [1, 0.7], color='black')

	#also draw the nice x-spine
	ax1.spines['left'].set_position('zero')
	ax1.spines['right'].set_color('none')
	ax1.spines['bottom'].set_position('zero')
	ax1.spines['top'].set_color('none')
	ax1.spines['left'].set_smart_bounds(True)
	ax1.spines['bottom'].set_smart_bounds(True)
	ax1.xaxis.set_ticks_position('bottom')

	plt.savefig(filename, format=imgfmt, dpi=dpi, bbox_inches='tight')


#OPTIMIZE ME: only parse once instead of every time

def summary(tab, html=False, outdir=None, prefix='', seqbank={}, tmcount={}, silent=False):
	out = ''
	if html == 1: out += '<pre>'
	#elif html == 2: out += '<table style="font-family: monospace, courier; font-size: 75%">'
	elif html == 2: out += '<table class="summtbl">'
	n = 0
	for l in tab.split('\n'):
		n += 1
		if not l.strip(): continue
		hit = l.split()[1]
		#if 'TC-DB' in hit: hit = hit.split('|')[3] + '-' + hit.split('|')[2]
		try: seq = seqbank[hit]
		except KeyError:
			seq = seqbank[hit] = subprocess.check_output(['blastdbcmd', '-db', 'tcdb', '-target_only','-entry', hit]).decode('utf-8')
		try: ntmss = tmcount[hit]
		except KeyError:
			tmss = parse_hmmtop(hmmtop(seq, outdir=outdir, silent=silent))
			ntmss = tmcount[hit] = len(tmss)

		tcid, acc = tuple(hit.split('-'))

		score = round(float((l.split()[11])))

		#query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

		fam = '%s.%s.%s' % tuple(tcid.split('.')[0:3])

		if html >= 2: 
			if n % 2: out += '<tr class="oddrow">'
			else: out += '<tr>'

		#COLUMN 1: Accession number/UNP
		if html >= 2: out += '<td>'

		if html: out += '<a href="http://www.tcdb.org/search/result.php?acc=%s">' % acc
		out += '%s' % acc
		if html: out += '</a>'
		if html >= 2: out += '</td>'

		#COLUMN 2: TMSs
		if html < 2: out += '\t'
		if html >= 2: out += '<td><a href="../blasts/%s_%s.top">' % (tcid, acc)
		out += '%d TMSs' % ntmss
		if html >= 2: out += '</a></td>'

		#COLUMN 3: TCID
		if html < 2: out += '\t'
		if html >= 2: out += '<td>'
		if html: out += '<a href="http://www.tcdb.org/search/result.php?tc=%s">' % fam
		out += '%s' % tcid
		if html: out += '</a>'
		if html >= 2: out += '</td>'

		#COLUMN 4: Name
		name = seq.split('\n')[0][1:]
		if html < 2:
			if len(name) > 37: name = name[:34] + '...'
		if html < 2: out += '\t'
		if html >= 2: out += '<td>'
		out += '%s' % name
		if html >= 2: out += '</td>'

		#COLUMN 5: Score
		if html < 2: out += '\t'
		if html >= 2: out += '<td>'
		if html: out += '<a href="#%s_%s">' % (prefix, acc)
		out += '%s' % score
		if html: out += '</a>'
		if html >= 2: out += '</td>'

		#COLUMN 6: evalue
		evalue = l.split()[10]
		if 'e' in evalue: evalue = evalue[evalue.find('e'):]
		if html < 2: out += '\t'
		if html >= 2: out += '<td>'
		out += '%s' % evalue
		if html >= 2: out += '</td>'

		if html >= 2: out += '</tr>'
		out += '\n'
	if html == 1: out += '</pre>'
	elif html == 2: out += '</table>'
	return out

def fmt_pairw(tab, pairw, html=0, prefix=''):
	pairw = pairw[pairw.find('\n>'):]
	pairw = pairw.split('Effective search space')[0]
	#re.sub('^>', '<hr/>\n>', pairw)

	final = ''
	redme = 0
	#records = pairw.split('\n>')[1:]
	#print(len(records), len(tab.split('\n')))
	#for x in zip(records, tab.split('\n')): 
	#	print(x[0].split()[0], x[1].split('\t')[1])

	name = ''
	names = []
	out = ''
	red = 0
	for l in pairw.split('\n'):
		if l.startswith('>'): 
			name = l.rstrip()

			if html: 
				#standard FASTA
				if '|' in name: acc = name.split('|')[2]
				#simplified TC-FASTA
				elif '-' in name: acc = name.split('-')[1].split(' ')[0]
				#worst case
				else: acc = name
				out += '<hr/><a name="%s_%s"></a>' % (prefix, acc )

		#only start after having found a name
		if name and l.startswith(' Score ='): 
			names.append(name)

		out += '\n'
		if html and red:
			out += '<span class="red">'
		out += l
		if html and red:
			red = 0
			out += '</span>'
		if l.startswith('Query'): red = 1
			
	#for l in zip(pairw.split('\n'), tab.split('\n')):
	#	print(l, l, l)
	#	print([l, l, l])
	#	if redme:
	#		redme -= 1
	#		final += '<span style="color:red">' + l[1] + '</span>'
	#	elif l[1].startswith('Query'): 
	#		redme = 1
	#		final += l[0]
	#	elif l[1].startswith('>'):
	#		final += '<a name=""'

	#pairw = pairw.replace('\n>', '<hr/>\n>')
	#return pairw
	return out

def til_warum(seq, outfile, title='Unnamed', dpi=100, html=2, outdir=None, clobber=False, seqbank={}, tmcount={}, silent=False, maxhits=50):

	tab, pairw = blast(seq, maxhits=maxhits)

	if not clobber and os.path.isfile(outfile):  pass
	else: plot_tab(tab, hmmtop(seq, outdir=outdir), outfile, dpi=dpi, overwrite=clobber, silent=silent)

	if not outfile.startswith('/'): relative = outfile[outfile.find('/')+1:]
	else: relative = outfile

	if outdir and not os.path.isdir(outdir): os.mkdir(outdir)

	out = ''

	#out += '<html><head><title>%s</title><link href="nice.css" type="text/css" rel="stylesheet"/></head><body>' % title
	outi = ''
	out += '<hr/>'

	out += summary(tab, html=html, outdir=outdir, prefix=title, tmcount=tmcount)
	out += '<hr/>'

	out += '<pre>'
	out += fmt_pairw(tab, pairw, html=2, prefix=title)
	#fmt_pairw(tab, pairw)
	out += '</pre>'
	#out += '</body></html>'
	

	return outi, out

if __name__ == '__main__': print('### UI NOT IMPLEMENTED, BUG KEVIN ###')
