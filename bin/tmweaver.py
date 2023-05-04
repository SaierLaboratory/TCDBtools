#!/usr/bin/env python

from __future__ import print_function, division
import os
import matplotlib.pyplot as plt

import argparse
import numpy as np

import sys
import subprocess
import re
import shlex

#TODO: generalize key classes and put them together
import libquod
import phoboshop

from Bio import SeqIO, Seq, AlignIO


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

def to_spans(indexlist):
    spans = []
    for i in sorted(indexlist):
        if not spans: spans.append([i, i])
        elif i == spans[-1][-1] + 1: spans[-1][-1] += 1
        else: spans.append([i, i])
    return spans

class Dictlike(dict):
    def get(self, index):
        if index in self: return self[index]
        elif index < min(self): return self[min(self)] + (index - min(self))
        elif index > max(self): return self[max(self)] + (index - max(self))
        else: raise Exception('impossible exception')

class Protein(phoboshop.Protein):
    pass
    vspans = []
    what = None
    spans = []

    hydro = None
    hmmtop = None

    def get_spans(self):
        #hmmtop = self.what.entities[1]
        #return hmmtop.spans
        return self.spans

    def pad(self):
        oldspans = self.spans
        offset = 0
        m = Dictlike()
        for resi, resn in enumerate(self.record.seq):
            if resn == '-': offset -= 1
            m[resi] = resi + offset
        
        self.spans = [
                list(np.clip([m.get(span[0]), m.get(span[-1])], 1, len(self.record)-1))
                for span in self.spans
                
                ]

    def plot_vspans(self, ax, style=0):
        color = TMCOLOR(style)
        alpha = 0.25

        vspans = []
        for span in self.spans:
            vspans.append(ax.axvspan(span[0], span[1], facecolor=color, alpha=alpha))
        self.vspans = vspans

    def render(self, style=0):
        self.hydro = self.entities['what'] = libquod.entities.Hydropathy.compute(self.record)
        self.hydro.set_edgecolor('r')
        #self.hmmtop = self.entities['hmmtop'] = libquod.entities.HMMTOP.compute(self.record)
        return [self.hydro]

    def hmmtop(self, tms=None):
        if tms is not None: 
            self.spans = tms
            return #libquod.entities.HMMTOP(spans=tms)
        fasta = '>untitled\n{}'.format(self.get_sequence())
        p = subprocess.Popen(['hmmtop', '-if=--', '-pi=spred', '-sf=FAS', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        try: out, err = p.communicate(input=fasta)
        except TypeError:
            p = subprocess.Popen(['hmmtop', '-if=--', '-pi=spred', '-sf=FAS', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            out, err = p.communicate(input=fasta.encode('utf-8'))
            out = out.decode('utf-8')


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
        hydro = self.hydro
        for x, y in zip(hydro.X, hydro.Y):
            if x == targetx: return 0 if np.isnan(y) else y
        return 0

class AlignmentBlock(object):
    seqlist = None
    def __init__(self, seqlist=None):
        self.seqlist = [] if seqlist is None else seqlist

    def __getitem__(self, index):
        return '\n'.join([seq[index] for seq in self.seqlist])

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
    occupancy = False
    def __init__(self, alignments):
        self.alignments = list(alignments)

    def get_header(self): return "MSA"

    def get_alignments(self): return self.alignments

    def get_sequence(self):
        for msa in self.alignments:
            seqlist = [str(record.seq) for record in msa]
            ablock = AlignmentBlock(seqlist)

        return ablock

    def hmmtop(self, tms=None):
        if tms is None: pass
        else: self.spans = tms

    def render(self, style=0):
        for msa in self.alignments:
            #self.what = quod.What('', style=style, nohmmtop=True)
            self.hydro = libquod.entities.BaseCurve(style=style)
            self.hydro.set_edgecolor('r')
            self.hydro.Y = avehas3.get_average_hydropathies(msa, window=19)
            self.hydro.X = np.arange(0., len(self.hydro.Y))

            

            self.entities['hydro'] = self.hydro
            entlist = []
            entlist = [self.hydro]

            if self.amphipathicity:
                amphi = libquod.entities.BaseCurve(style='g')
                amphi.Y = avehas3.get_average_amphipathicities(msa, window=19)
                amphi.X = np.arange(0, len(amphi.Y))
                entlist.append(amphi)

            entlist.append(avehas3.TMcenter(avehas3.get_tmcenters(msa), ymin=-3, ymax=-2.4))

            simil = avehas3.get_similarities(msa)
            window = 10
            ymin = -3
            ymax = -1.5
            scores = np.array([np.nanmean(simil[i:i+window]) for i in range(0, len(simil)-window)])
            normscores = (scores - min(scores)) / (max(scores) - min(scores))
            scaledscores = normscores * (ymax - ymin)
            shiftedscores = scaledscores + ymin
            similcurve = libquod.entities.BaseCurve(style='gray')
            similcurve.Y = shiftedscores
            similcurve.X = np.arange(0, len(shiftedscores)) + window//2
            similcurve.set_edgecolor('gray')
            entlist.append(similcurve)

            if self.occupancy:
                occcurve = libquod.entities.BaseCurve('', style='y')
                occcurve.Y = avehas3.get_occupancies(msa)
                occcurve.X = np.arange(len(occcurve.Y))
                entlist.append(occcurve)
            

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
        self.hydro = libquod.entities.Hydropathy.compute(seq=self.fullfasta, fragment=self.fragfasta, style=style)
        self.hydro.set_edgecolor('r')
        self.entities['hydro'] = self.hydro
        entlist = []
        entlist = [self.hydro]

        return entlist

    def hmmtop(self, tms=None):
        if tms is not None: 
            self.spans = tms
            return
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

    def add_protein(self, protein, tms=None):
        self.proteins.append(protein)
        self.tmss.append(protein.hmmtop(tms))

    def add_alignment(self, alignment, tms):
        self.proteins.append(alignment)
        alignment.hmmtop(tms)

    def run(self):
        print('Use [?] to get help on phoboshop shortcuts', file=sys.stderr)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        for i, p in enumerate(self.proteins):
            for e in p.render(style=i):
                e.plot(self.ax)
            p.plot_vspans(self.ax, style=i)

        self.update_title()

        self.fig.canvas.mpl_connect('button_press_event', self.onmousedown)
        self.fig.canvas.mpl_connect('button_release_event', self.onmouseup)
        self.fig.canvas.mpl_connect('motion_notify_event', self.onmousemove)
        self.fig.canvas.mpl_connect('key_press_event', self.onkeydown)

        self.testline, = self.ax.plot([], [])
        self.ax.axhline(0, color='k', lw=0.5)
        self.ax.set_ylim([-3, 3])
        plt.show()

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
            vspan = self.ax.axvspan(event.xdata, event.xdata+1, fc=TMCOLOR(self.pid), alpha=0.25)
            if self.pid is None: 
                self.proteins[0].vspans.append(vspan)
                self.target = (0, len(self.proteins[0].vspans)-1)
            else: 
                self.proteins[self.pid].vspans.append(vspan)
                self.target = (self.pid, len(self.proteins[self.pid].vspans)-1)
            self.ax.figure.canvas.draw()

        elif self.mode.startswith('mark'):
            self.mode += 'drag'
            target = []
            for i, p in enumerate(self.proteins):
                x = int(round(event.xdata))
                y = p.get_ydro(x)
                marker = self.ax.scatter([x], [y], c='k')
                target.append(marker)
                try: resn = p.get_sequence()[x]
                except IndexError: resn = ''
                label = self.ax.text(x, y, resn)
                target.append(label)
            self.target = target
            
            self.fig.canvas.draw()

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
                print('#{}{:d} in {} (= {:0.2f}kcal/mol)'.format(label.get_text(), int(round(marker.get_offsets()[0,0])), self.proteins[n].get_header(), marker.get_offsets()[0,1]), file=sys.stderr)
                marker.remove()
                label.remove()
                n += 1
            self.target = None
            self.fig.canvas.draw()

        elif 'delete' in self.mode:
            inspan = self._is_in_span(event.xdata, pid=None)
            if not inspan: return

            vspan = self.proteins[inspan[0]].vspans.pop(inspan[1])
            self.proteins[inspan[0]].spans.pop(inspan[1])
            vspan.remove()
            self.fig.canvas.draw()
            
    def _cleave(self, string, suffix):
        if suffix in string:
            return string[:string.find(suffix)]
        else: return string

    def onmousemove(self, event): 
        if 'drag' not in self.mode: return
        elif event.inaxes != self.ax: return

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
            self.fig.canvas.draw()

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

        elif event.key == 'q':
            sys.exit(0)

        elif event.key == 'r':
            self.rotate_format()

        elif event.key == 'w':
            self.write_tmss()

        elif event.key == 'W':
            self.write_tmss_relative()

        elif event.key == 'escape':
            self.mode = 'normal'

        self.update_title()

    def rotate_format(self):
        cycle = {'fasta':'hmmtop', 'hmmtop':'multiquod', 'multiquod':'ranges', 'ranges':'tsv', 'tsv':'fasta'}
        print('Switched outfmt from {} to {}'.format(self.outfmt, cycle[self.outfmt]))
        self.outfmt = cycle[self.outfmt]

    def write_protein(self, protein):
        out = ''
        if self.outfmt == 'fasta':
            s = ''
            for i, span in enumerate(sorted(protein.get_spans())):
                s += '>{}_TMS{}\n'.format(protein.record.id, i+1)
                s += str(protein.get_sequence()[span[0]-1:span[1]])
                s += '\n'
            out += s + '\n'
        elif self.outfmt == 'hmmtop':
            s = '>HP: {} {}  UNK  {}  '.format(len(protein.get_sequence()), protein.record.id, len(protein.get_spans()))
            for span in sorted(protein.get_spans()):
                s += '{}  {}  '.format(span[0], span[1])
            out += '{}\n'.format(s)
        elif self.outfmt == 'multiquod':
            s = '#add TMSs for {}\n'.format(protein.record.id)
            s += 'tms add SUBPLOT COLOR '
            for span in sorted(protein.get_spans()):
                s += ' {} {}'.format(span[0], span[1])
            out += '{}\n'.format(s)
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
        else: raise NotImplementedError('Unknown outfmt: {}'.format(self.outfmt))
        return out

    def write_tmss(self, pad=False):
        out = ''
        for protein in self.proteins:
            if type(protein) is Protein:
                out += self.write_protein(protein)
            elif type(protein) is Alignment:
                for alignment in protein.alignments:
                    for record in alignment:
                        p = Protein(record)
                        for e in protein.entities:
                            p.spans = protein.spans
                            if pad: p.pad()
                        out += self.write_protein(p)
            else: raise TypeError('Unrecognized type')
            out = out.strip()

            if self.append: 
                with open(self.outfile, 'a') as f: 
                    f.write(out + '\n')
                    if self.outfile != "/dev/stdout":
                        print("Appended {} output to {}".format(self.outfmt, self.outfile), file=sys.stderr)
            else:
                with open(self.outfile, 'w') as f: 
                    f.write(out + '\n')
                    if self.outfile != "/dev/stdout":
                        print("Wrote {} output to {}".format(self.outfmt, self.outfile), file=sys.stderr)

    def write_tmss_relative(self):
        return self.write_tmss(pad=True)

        #out = ''
        #for protein in self.proteins:
        #    seqspanlists = []
        #    seqrecords = []
        #    if isinstance(protein, Alignment):
        #        allspans = set()
        #        for span in protein.get_spans():
        #            allspans = allspans.union(set([i for i in range(span[0], span[1]+1)]))

        #        for msa in protein.alignments:
        #            for record in msa:
        #                seqrecords.append(record)
        #                seqresi = 0
        #                seqindices = set()
        #                for resi, resn in enumerate(record.seq):
        #                    if resi in allspans: seqindices.add(seqresi)

        #                    if resn in 'ACDEFGHIKLMNPQRSTVWY': seqresi += 1
        #                seqspanlists.append(to_spans(seqindices))

        #    else: #assume it's a sequence-like
        #        allspans = set()
        #        for span in protein.get_spans():
        #            allspans = allspans.union(set([i for i in range(span[0], span[1]+1)]))

        #        seqrecords.append(protein.record)
        #        seqresi = 0
        #        seqindices = set()
        #        for resi, resn in enumerate(protein.record.seq):
        #            if resi in allspans: seqindices.add(seqresi)
        #            if resn in 'ACDEFGHIKLMNPQRSTVWY': seqresi += 1

        #        seqspanlists.append(to_spans(seqindices))

        #    if self.outfmt == 'fasta':
        #        for record, spans in zip(seqrecords, seqspanlists):
        #            s = ''
        #            for spani, span in enumerate(sorted(spans)):
        #                s += '>{}_TMS{}\n'.format(record.id, spani+1)
        #                s += record.seq[span[0]-1:span[1]] + '\n'
        #            out += '{}\n'.format(s)

        #    elif self.outfmt == 'hmmtop':
        #        for record, spans in zip(seqrecords, seqspanlists):
        #            s = '>HA: {} {}  UNK  {}  '.format(len(record.seq), record.id, len(spans))
        #            for span in sorted(spans):
        #                s += '{}  {}  '.format(span[0], span[-1])
        #            out += '{}\n'.format(s)

        #    elif self.outfmt == 'multiquod':
        #        for record, spans in zip(seqrecords, seqspanlists):
        #            s = '#add TMSs for {}\n'.format(record.id)
        #            s += 'tms add SUBPLOT COLOR '
        #            for span in sorted(spans):
        #                s += ' {} {}'.format(span[0], span[-1])
        #            out += '{}\n'.format(s)

        #    elif self.outfmt == 'ranges':
        #        for record, spans in zip(seqrecords, seqspanlists):
        #            s = '{}\t'.format(record.id)
        #            for span in sorted(spans):
        #                s += '{}-{},'.format(span[0], span[-1])
        #            if s.endswith(','): s = s[:-1]
        #            out += s + '\n'

        #    elif self.outfmt == 'tsv':
        #        for record, spans in zip(seqrecords, seqspanlists):
        #            s = '{}\t'.format(record.id)
        #            for span in spans:
        #                for x in span: s += '{}\t'.format(x)
        #            out += s.strip() + '\n'
        #    else: raise NotImplementedError('Unknown outfmt: {}'.format(self.outfmt))

        #out = out.strip()

        #if self.append: 
        #    with open(self.outfile, 'a') as f: f.write(out + '\n')
        #else:
        #    with open(self.outfile, 'w') as f: f.write(out + '\n')

    def update_title(self):
        title = 'Insert title here'
        if self.mode == 'normal': title = 'NORMAL mode'
        elif self.mode.startswith('edit'): 
            if self.pid is None: title = 'EDIT mode'
            else: title = 'EDIT {} mode'.format(self.pid+1)
        elif self.mode.startswith('delete'): title = 'DELETE mode'
        elif self.mode.startswith('mark'): title = 'MARK mode'

        plt.get_current_fig_manager().set_window_title(title)

    def print_help(self):
        print('''
CONTROLS
========

1-9 - Enter edit mode for a specific sequence
a - Enter edit mode
d - Enter deletion mode
i - Enter edit mode
m - Enter mark mode
r - Switch output format
W - Write cut sequence(s) to disk/stdout with relative indices. TODO: extend to FASTAs, run the selected pipeto program if defined (useful for TMSOC, presumably)
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

    parser.add_argument('infile', nargs='+', help='Sequence files to read')
    parser.add_argument('-o', metavar='OUTFILE', default='/dev/stdout', help='Where to save the new indices. (default: print to stdout)')
    parser.add_argument('-m', '--mode', default='multiquod', help='Format to print in indices in. (fasta, hmmtop, \033[1mmultiquod\033[0m, ranges, tsv)')
    parser.add_argument('--amphipathicity', action='store_true', help='Also plot amphipathicity (currently only implemented for MSAs)')
    parser.add_argument('--occupancy', action='store_true', help='Also plot occupancy (currently only implemented for MSAs)')
    parser.add_argument('-a', action='store_true', help='Append to file instead of overwriting')
    parser.add_argument('--frag', action='store_true', help='Open sequences as fragments (frag1 full1 frag2 full2... fragN fullN)')
    parser.add_argument('-t', '--tms', nargs='+', help='Space-separated dash-delimited ranges representing precomputed TMSs. Bypasses HMMTOP if set.')
    #TODO: decide how to handle overlapping "TMSs"
    args = parser.parse_args()

    tmweaver = TMWeaver()
    tmweaver.outfile = args.o
    tmweaver.append = args.a

    modes = ('fasta', 'hmmtop', 'multiquod', 'ranges', 'tsv')
    if args.mode not in modes: 
        parser.print_help()
        exit(1)
    else: tmweaver.outfmt = args.mode

    if args.tms is not None:
        tms = [[int(x) for x in span.split('-')] for span in args.tms]
    else: tms = None

    if args.frag:
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
            tmweaver.add_protein(p, tms)
    else:
        for infn in args.infile:
            f = open(infn)
            filetype = detect_filetype(f)
            f.seek(0)
            if filetype == 'fasta':
                for record in SeqIO.parse(f, 'fasta'):
                    tmweaver.add_protein(Protein(record), tms)
                f.close()
            elif filetype == 'clustal':
                if avehas3 is None: raise ImportError('Could not find AveHAS3')
                alignments = avehas3.AlignIO.parse(f, 'clustal')
                aln = Alignment(alignments=alignments)
                aln.amphipathicity = args.amphipathicity
                aln.occupancy = args.occupancy
                aln.header = infn
                tmweaver.add_alignment(aln, tms)

    tmweaver.run()
