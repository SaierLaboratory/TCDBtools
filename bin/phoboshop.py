#!/usr/bin/env python

from __future__ import print_function, division
import argparse

import os
import matplotlib.pyplot as plt
import libquod
import numpy as np

import sys
import subprocess
import re
import shlex

from Bio import SeqIO, Seq


class Protein(object):
    fasta = ''
    record = None
    entities = {}
    def __init__(self, record):
        self.record = record

    def __hash__(self):
        return hash((type(self), self.record.id, self.record.seq)) #please don't use MutableSeq

    def get_header(self): return self.record.id

    def get_entropy(self, a=None, b=None):
        seq = self.get_sequence(a, b)

        #seq = seq[int(np.floor(a)):int(np.ceil(b))]
        total = 0
        counts = {}
        for resn in seq:
            try: counts[resn] += 1
            except KeyError: counts[resn] = 1
            total += 1
        entropy = 0
        for k in counts: entropy -= counts[k]/total*np.log2(counts[k]/total)
        return entropy

    def get_sequence(self, a=None, b=None):
        return self.record.seq[a if a is None else max(0, int(np.floor(a))):b if b is None else min(len(self.seq)-1, int(np.floor(b)))]

    def get_header(self):
        return '>{}'.format(self.record.id)

    def render(self, style=0):
        self.entities['hydro'] = libquod.entities.Hydropathy.compute(self.record)
        self.entities['hmmtop'] = libquod.entities.HMMTOP.compute(self.record)

        entlist = []
        for k in self.entities: entlist.append(self.entities[k])
        return entlist

class Phoboshop(object):

    proteins = []
    mode = 'normal'
    fig = None
    plot = None

    outfile = '/dev/stdout'

    pipeto = None

    testline = None
    selections = {}
    queue = {}

    def __init__(self):
        pass

    def add_protein(self, protein):
        self.proteins.append(protein)

    def run(self):
        print('Use [?] to get help on phoboshop shortcuts', file=sys.stderr)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.entities = {}
        for i, p in enumerate(self.proteins): 
            self.entities[p] = []
            for e in p.render(style=i): 
                self.entities[p].append(e)
                e.plot(self.ax)

        self.update_title()
        self.fig.canvas.mpl_connect('button_press_event', self.onmousedown)
        self.fig.canvas.mpl_connect('button_release_event', self.onmouseup)
        self.fig.canvas.mpl_connect('motion_notify_event', self.onmousemove)
        self.fig.canvas.mpl_connect('key_press_event', self.onkeydown)
        self.testline, = self.ax.plot([], [])

        plt.show()

    def onmousedown(self, event):
        if not event.inaxes: return

        if self.mode.startswith('cut'):
            self.mode = 'cutdrag'
            if 'cut' in self.selections and self.selections['cut'] is not None:
                vspan = self.selections['cut']
            else: 
                vspan = self.ax.axvspan(event.xdata, event.xdata, fc='red', alpha=0.4, hatch='/')
                self.selections['cut'] = vspan
            xy = vspan.get_xy()
            xy[:,0] = event.xdata
            xy[2:4,0] = event.xdata + 1
            vspan.set_xy(xy)
            vspan.figure.canvas.draw()

        elif self.mode == 'delete':
            if 'cut' in self.queue:
                popme = []
                for i, vspan in enumerate(self.queue['cut']):
                    xy = vspan.get_xy()
                    xmin = min(xy[0,0], xy[2,0])
                    xmax = max(xy[0,0], xy[2,0])
                    if xmin <= event.xdata <= xmax: 
                        vspan.remove()
                        popme.append(i)
                        self.fig.canvas.draw()
                for i in popme[::-1]: self.queue['cut'].pop(i)
            if 'cut' in self.selections and self.selections['cut'] is not None:
                vspan = self.selections['cut']
                xy = vspan.get_xy()
                xmin = xy[0,0]
                xmax = xy[2,0]
                if xmin <= event.xdata <= xmax: 
                    vspan.remove()
                    del vspan
                    self.fig.canvas.draw()
            if 'divide' in self.queue:
                popme = []
                for i, vline in enumerate(self.queue['divide']):
                    x = vline.get_xdata()
                    if x[0]-2 <= event.xdata <= x[0]+2:
                        vline.remove()
                        popme.append(i)
                if popme:
                    for i in popme[::-1]: self.queue['divide'].pop(i)
                    self.fig.canvas.draw()
            if 'divide' in self.selections and self.selections['divide'] is not None:
                vline = self.selections['divide']
                x = vline.get_xdata()
                if x[0]-2 <= event.xdata <= x[0]+2:
                    vline.remove()
                    del vline
                    self.fig.canvas.draw()


        elif self.mode == 'divide':
            self.mode = 'dividedrag'
            if 'divide' in self.selections and self.selections['divide'] is not None:
                vline = self.selections['divide']
            else:
                vline = self.selections['divide'] = self.ax.axvline(x=int(event.xdata), color='k')
            x = vline.get_xdata()
            x[0] = event.xdata
            x[1] = event.xdata
            vline.set_xdata(x)
            
            vline.figure.canvas.draw()

        elif self.mode == 'entropy':
            vspan = self.ax.axvspan(event.xdata, event.xdata+1, fc='g', alpha=0.2, hatch='\\')
            self.mode = 'entropydrag'
            self.selections['entropy'] = vspan
            vspan.figure.canvas.draw()

        elif self.mode == 'move':
            if 'divide' in self.queue:
                for vline in self.queue['divide']:
                    x = vline.get_xdata()
                    if x[0]-2 <= event.xdata <= x[0]+2:
                        self.mode = 'movedrag'
                        self.selections['move'] = (vline, 'divide', event.xdata)
                        return
            if 'cut' in self.queue:
                for vspan in self.queue['cut']:
                    xy = vspan.get_xy()
                    
                    xmin = min(xy[0,0], xy[2,0])
                    xmax = max(xy[0,0], xy[2,0])
                    #edges
                    if xmax-2 <= event.xdata <= xmax+2:
                        self.mode = 'moveresizedrag'
                        self.selections['resize'] = (vspan, 'cut', 'r')
                        return
                    elif xmin-2 <= event.xdata <= xmin+2:
                        self.mode = 'moveresizedrag'
                        self.selections['resize'] = (vspan, 'cut', 'l')
                        return
                    elif xmin <= event.xdata <= xmax:
                        self.mode = 'movedrag'
                        self.selections['move'] = (vspan, 'cut', event.xdata)
                        return
            
    def onmouseup(self, event):
        #if self.mode.startswith('cut'): self.mode = 'cut'
        if self.mode == 'cutdrag':
            try: self.queue['cut']
            except KeyError: self.queue['cut'] = []

            vspan = self.selections['cut']
            xy = vspan.get_xy()
            xmin = min(xy[:,0])
            xmax = max(xy[:,0])
            xy[:,0] = xmin
            xy[2:4,0] = xmax

            self.queue['cut'].append(self.selections['cut'])
            self.selections['cut'] = None
            self.mode = 'cut'
        elif self.mode == 'dividedrag':
            try: self.queue['divide']
            except KeyError: self.queue['divide'] = []

            vline = self.selections['divide']
            self.queue['divide'].append(self.selections['divide'])
            self.selections['divide'] = None
            self.mode = 'divide'

        elif self.mode == 'entropydrag':
            vspan = self.selections['entropy']

            xy = vspan.get_xy()
            xmin = int(np.floor(min(xy[:,0])))
            xmax = int(np.ceil(max(xy[:,0])))
            for p in self.proteins:
                e = p.get_entropy(xmin, xmax)
                length = len(p.get_sequence(xmin, xmax))
                if length == 0:
                    print('{}: {:0.3f} bits/aa({}-{}: {} aa, max entropy: {:0.3f} bits/aa)'.format(
                        p.record.id,
                        e, 
                        xmin, xmax, 
                        length,
                        min(-np.log2(1/20), 0)
                    ))
                else:
                    print('{}: {:0.3f} bits/aa({}-{}: {} aa, max entropy: {:0.3f} bits/aa)'.format(
                        p.record.id,
                        e, 
                        xmin, xmax, 
                        length,
                        min(-np.log2(1/20), -np.log2(1/(length)))
                    ))
                #print(p.get_sequence(xmin, xmax))
                

            vspan.remove()
            del vspan
            self.selections['entropy'] = None
            self.mode = 'entropy'
            self.fig.canvas.draw()

        elif self.mode == 'movedrag':
            self.selections['move'] = None
            self.mode = 'move'
        elif self.mode == 'moveresizedrag':

            vspan, side = self.selections['resize']
            xy = vspan.get_xy()
            xmin = min(xy[:,0])
            xmax = max(xy[:,0])
            xy[:,0] = xmin
            xy[2:4,0] = xmax

            self.selections['resize'] = None
            self.mode = 'move'
        #if self.mode.startswith('cut'): self.mode = 'cut'


    def onmousemove(self, event):
        if self.mode == 'cutdrag':
            if event.xdata is None: return
            if 'cut' in self.selections:
                vspan = self.selections['cut']
                xy = vspan.get_xy()
                xy[2:4,0] = event.xdata
                vspan.set_xy(xy)
                vspan.figure.canvas.draw()

        elif self.mode == 'dividedrag':
            if event.xdata is None: return
            if 'divide' in self.selections:
                vline = self.selections['divide']
                x = vline.get_xdata()
                x[0] = int(event.xdata)
                x[1] = int(event.xdata)
                vline.set_xdata(x)
                self.fig.canvas.draw()

        elif self.mode == 'entropydrag':
            if event.xdata is None: return
            vspan = self.selections['entropy']
            xy = vspan.get_xy()
            xy[2:4,0] = event.xdata
            vspan.set_xy(xy)
            vspan.figure.canvas.draw()

        elif self.mode == 'movedrag':
            if event.xdata is None: return
            vspan, t, orig = self.selections['move']
            if t == 'cut':
                xy = vspan.get_xy()
                xy[:,0] += (event.xdata - orig)
                vspan.set_xy(xy)
                vspan.figure.canvas.draw()
                self.selections['move'] = (vspan, t, event.xdata)
            elif t == 'divide':
                x = vspan.get_xdata()
                x[0] += event.xdata - orig
                x[1] += event.xdata - orig
                vspan.set_xdata(x)
                self.fig.canvas.draw()
                self.selections['move'] = (vspan, t, event.xdata)

        elif self.mode == 'moveresizedrag':
            if event.xdata is None: return
            vspan, side = self.selections['resize']
            xy = vspan.get_xy()
            if side == 'l':
                xy[:2,0] = event.xdata
                xy[4:,0] = event.xdata
            elif side == 'r':
                xy[2:4,0] = event.xdata
            vspan.set_xy(xy)
            vspan.figure.canvas.draw()
        else: return


    def _merge(self, spanlist):
        def _overlap(span1, span2):
            # |---|
            #  |-|
            if (span1[0] <= span2[0]) and (span2[1] <= span1[1]): return True
            #  |-|
            # |---|
            elif (span2[0] <= span1[0]) and (span1[1] <= span2[1]): return True
            #  |---|
            # |---|
            elif (span2[0] <= span1[0] <= span2[1]): return True
            elif (span1[0] <= span2[1] <= span1[1]): return True
            # |---|
            #  |---|
            elif (span1[0] <= span2[0] <= span1[1]): return True
            elif (span2[0] <= span1[1] <= span2[1]): return True
            else: return False
        def _check_overlaps(spanlist):
            for i, span1 in enumerate(spanlist):
                for span2 in spanlist[i+1:]:
                    if _overlap(span1, span2): return True
            return False

        spanlist.sort()
        while _check_overlaps(spanlist):
            popme = []
            for i, span1 in enumerate(spanlist):
                for j, span2 in enumerate(spanlist):
                    if j <= i: continue
                    if _overlap(span1, span2): 
                        span1[1] = max(span1[1], span2[1])
                        popme.append(j)
            #for j in popme[::-1]: spanlist.pop(j)
            spanlist.pop(popme[-1])
        return spanlist
                    
        

    def write_cuts(self):
        maxlen = max([len(p.record) for p in self.proteins])
        indices = set(range(maxlen))
        if 'divide' in self.queue:
            divpoints = list(reversed(sorted([int(vline.get_xdata()[0]) for vline in self.queue['divide']])))
            indexbins = {}
            binid = 0
            for i in range(maxlen):
                if divpoints and i >= divpoints[-1]:
                    divpoints.pop()
                    binid += 1
                indexbins[i] = binid
                
        if 'cut' in self.queue:
            for vspan in self.queue['cut']:
                xy = vspan.get_xy()
                xmin = int(round(min(xy[:,0])))
                xmax = int(round(max(xy[:,0])))

                for x in range(xmin, xmax+1):
                    if x in indices: indices.remove(x)

        
        s = ''
        for p in self.proteins:
            if 'divide' in self.queue:
                seqlist = {}
                for resi in indexbins: seqlist[indexbins[resi]] = ''

                for resi in indices:
                    if 0 <= resi < len(p.record):
                        seqlist[indexbins[resi]] += p.record.seq[resi]
                for seqid in seqlist:
                    s += format(SeqIO.SeqRecord(id='{}_part{}'.format(p.record.id, seqid), name='', description='', seq=Seq.Seq(seqlist[seqid])), 'fasta') if seqlist[seqid] else ''
            else:
                s += p.get_header() + '\n'
                seq = p.get_sequence()
                s += ''.join([seq[_] for _ in indices if (0 <= _ < len(seq))])
            s += '\n'
        s += '\n'
        if self.pipeto:
            p = subprocess.Popen(shlex.split(self.pipeto), stdin=subprocess.PIPE)
            out, err = p.communicate(input=s)
        with open(self.outfile, 'w') as f:
            f.write(str(s))
                

    def print_help(self):
        print('''
CONTROLS
========

c - Enter cutting mode. Click and drag to mark unwanted segments
d - Enter selection-deleting mode. As a double negative is a positive, clicking on segments marked for deletion marks them as wanted once again
e - Enter entropy mode. Click and drag to probe for local sequence complexity
m - Enter move/resize mode. Click and drag cores of selections to move them or edges to resize them
v - Enter diVision mode. Click to place dividers
w - Write cut sequence(s) to disk/stdout and run the selected pipeto program if defined
? - Pull up this help page
ESC - Enter normal mode, which does absolutely nothing
    ''', file=sys.stderr)

    def onkeydown(self, event):
        #prevent mode-switching while dragging stuff around
        if 'drag' in self.mode: return

        if event.key == 'c':
            self.mode = 'cut'
        elif event.key == 'd':
            self.mode = 'delete'
        elif event.key == 'e':
            self.mode = 'entropy'
        elif event.key == 'escape':
            self.mode = 'normal'

        elif event.key == 'm':
            self.mode = 'move'

        elif event.key == 'v':
            self.mode = 'divide'

        elif event.key == 'w':
            self.write_cuts()

        elif event.key == '?':
            self.print_help()

        #elif event.key == 'm':
        #    if event.inaxes != self.plot.ax: return
        #    print(event.xdata, event.ydata)
        #    xs = list(self.testline.get_xdata())
        #    xs.append(event.xdata)
        #    self.testline.set_xdata(xs)
        #    ys = list(self.testline.get_ydata())
        #    ys.append(event.ydata)
        #    self.testline.set_ydata(ys)
        #    self.testline.figure.canvas.draw()
        #    #self.plot.ax.

        self.update_title()

    def update_title(self):
        if self.mode == 'normal': plt.get_current_fig_manager().set_window_title('Normal mode')
        elif self.mode.startswith('cut'): plt.get_current_fig_manager().set_window_title('SELECT-FOR-DELETION mode')
        elif self.mode.startswith('divide'): plt.get_current_fig_manager().set_window_title('DIVIDE mode')
        elif self.mode.startswith('move'): plt.get_current_fig_manager().set_window_title('MOVE-SELECTION mode')
        elif self.mode.startswith('entropy'): plt.get_current_fig_manager().set_window_title('ENTROPY mode')
        elif self.mode == 'delete': plt.get_current_fig_manager().set_window_title('DESELECT mode')

def split_fasta(f):
    firstline = True
    sequences = []
    for l in f:
        if firstline and not l.startswith('>'): sequences.append('>sequence\n')
        firstline = False
        if l.startswith('>'):
            sequences.append(l)
        else:
            sequences[-1] += l
    return sequences

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('infile', nargs='+', help='Sequence files to read')
    parser.add_argument('-o', metavar='OUTFILE', default='/dev/stdout', help="Where to save the cut sequence(s). (default: print to stdout)")
    parser.add_argument('-p', metavar='PROGNAME | "PROGNAME ARG1 ARG2..."', help="Where to pipe the output (in addition to saving the cut sequence(s). Specify `xclip' to send stuff to the clipboard on Linux and `pbcopy' to send stuff to the clipboard on Mac. If the target program requires multiple arguments, enclose PROGNAME in quotes")

    args = parser.parse_args()


    phoboshop = Phoboshop()
    phoboshop.outfile = args.o

    if args.p: phoboshop.pipeto = args.p

    for infn in args.infile:
        for record in SeqIO.parse(infn, 'fasta'):
            p = Protein(record)
            phoboshop.add_protein(p)
    phoboshop.run()
