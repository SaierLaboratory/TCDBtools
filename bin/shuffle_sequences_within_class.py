#!/usr/bin/env python
from __future__ import print_function

import argparse
import random
import subprocess
import re
import pathlib
import Bio.SeqIO, Bio.Seq, Bio.AlignIO, Bio.Align
import numpy as np
import pathlib
import sys

from kevtools_common.parsing import parse_ranges

def hmmtop(seq):
    ''' Run HMMTOP on SeqRecord objects '''

    mapping = {}
    offset = 1 #HMMTOP is 1-indexed
    for trueresi, resn in enumerate(seq):
        if resn == '-': offset -= 1
        #issue: collapsed gap indices?
        #favor N:
        if trueresi + offset not in mapping: mapping[trueresi + offset] = trueresi 
        ##favor C:
        #mapping[trueresi + offset] = trueresi 
        
    
    topo = None
    indices = []

    p = subprocess.Popen(['hmmtop', '-if=--', '-is=pseudo', '-pi=spred', '-sf=FAS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = p.communicate(input=format(seq, 'fasta').encode('utf-8'))

    out = out.decode('utf-8')

    if out.strip():
        relevant = re.findall('(?:IN|OUT)\s+(?:[0-9]+\s*)+$', out)[0]
        if relevant.startswith('IN'): topo = 1
        elif relevant.startswith('OUT'): topo = 0
        else: topo = None
    else: 
        topo = 1 #because most proteins end up in the cytosol
        relevant = ''
    x = relevant.split()
    #indices = [[int(x[i])-1, int(x[i+1])-1] for i in range(2, len(x), 2)]
    indices = [[mapping[int(x[i])], mapping[int(x[i+1])]] for i in range(2, len(x), 2)]

    return topo, indices

def gen_states(topo, indices, seqlen):
    ''' Fill out per-position states given initial topology, TMS indices, and sequence length

    topo: Initial state represented as 0 for non-cytosolic environments (i.e. `O`utside) or 1 for cytosolic environments (i.e. `I`nside).
    indices: List of lists containing the positions of TMS boundary residues.
    seqlen: Total length of sequence. This is used to pad out the states so they match the original sequence in length.
    '''

    if topo == 0:
        start = 'O'
        inverse = 'I'
    elif topo == 1:
        start = 'I'
        inverse = 'O'
    else: raise ValueError

    states = [start for i in range(seqlen)]

    for i, tms in enumerate(indices):
        for j in range(tms[0], tms[1]+1): 
            states[j] = 'H'

    for i, tms in enumerate(indices):
        for j in range(tms[0]-1, 0, -1): 
            if states[j] == 'H': break
            elif (i % 2): #before an even (note 1 vs. 0-indexing) TMS: inverse
                states[j] = inverse
            elif not (i % 2): #before an odd (note 1 vs. 0-indexing) TMS:normal 
                states[j] = start
        for j in range(tms[1]+1, seqlen, 1):
            if states[j] == 'H': break
            elif (i % 2): #before an even (note 1 vs. 0-indexing) TMS: inverse
                states[j] = start 
            elif not (i % 2): #before an odd (note 1 vs. 0-indexing) TMS:normal 
                states[j] = inverse

    for i, tms in enumerate(indices):
        for j in range(tms[1]+1, min(seqlen, tms[1]+1+15), 1):
            if states[j] == 'H': break
            else: states[j] = states[j].lower()

        for j in range(tms[0]-1, max(0, tms[0]-1-15), -1): 
            if states[j] == 'H': break
            else: states[j] = states[j].lower()

    return states


def get_shuffle(template, mergedict=None):
    ''' Deprecated legacy function '''
    shuffleme = {}
    for k in template:
        if not mergedict:
            shuffleme[k] = template[k][:]
        else:
            if k in mergedict:
                shuffleme[k] = shuffleme[mergedict[k]]
            else:
                shuffleme[k] = template[k][:]
    return shuffleme


def get_pieces(seq, mergedict=None, topo=None, indices=None):
    ''' Get a shuffleable representation of a sequence and its state sequence. 

    seq: The original sequence. For sequence-agnostic cases, a range() with the appropriate length suffices. 
    mergedict: Dictionary mapping equivalent states together.
    topo: Initial state as represented by an integer. 0=outer/extracellular/luminal, 1=inner/cytoplasmic.
    indices: List of lists containing TMS boundaries.
    '''
    if indices is None: topo, indices = hmmtop(seq)
    if topo is None: topo = 1
    states = gen_states(topo, indices, len(seq))
    pieces = {}
    for k in 'HIiOo': pieces[k] = ''
    for resn, state in zip(seq, states):
        pieces[state] += resn

    template = {}
    for k in 'HIiOo':
        if mergedict: 
            if k in mergedict:
                try: template[mergedict[k]].extend(list(pieces[k]))
                except KeyError: template[mergedict[k]] = list(pieces[k])
                template[k] = template[mergedict[k]]
            else: template[k] = list(pieces[k])
        else: template[k] = list(pieces[k])
        #template[k] = list(pieces[k])
    return template, states

def old_shuffle_seq(seq, count=1, prefix='shuffled', mergedict=None, topo=None, indices=None):
    ''' Deprecated legacy function '''
    sequences = []
    template, states = get_pieces(seq, mergedict=mergedict, topo=topo, indices=indices)
    for i in range(count):
        shuffleme = get_shuffle(template, mergedict=mergedict)
        for k in shuffleme:
            if k in mergedict: continue
            random.shuffle(shuffleme[k])
        newseq = ''

        #for k in shuffleme: print(k, len(template[k]), len(shuffleme[k]), states.count(k))

        for state in states:
            newseq += shuffleme[state].pop(0)
        sequences.append('>{}_{:016X}\n{}\n'.format(prefix, abs(hash(newseq)), newseq))
        del shuffleme
    return sequences

def shuffle_seq(seq, count=1, prefix='shuffled', mergedict=None, topo=None, indices=None):
    ''' Wrapper for running a one-off shuffler 
    seq: SeqRecord to shuffle 
    count: '''
    template, states = get_pieces(seq, mergedict=mergedict, topo=topo, indices=indices)
    shuffler = StateShuffler(mergedict, states)
    return [Bio.SeqIO.SeqRecord(
            id='{}_{:016X}'.format(prefix, abs(hash(''.join(newseq)))), 
            name='', description='', 
            seq=Bio.Seq.Seq(''.join(newseq)),
        ) for newseq in shuffler.shuffled(seq.seq, count=count)]

class StateShuffler(object):
    ''' Shuffler object that can generate many within-state shuffles given a sequence of states. '''
    def __init__(self, mergedict, states):
        ''' Constructor for StateShuffler

        mergedict: Dictionary with equivalent state mappings. Residues belonging to states marked as identical in mergedict are eligible for shuffling with each other. 
        states: List of per-residue states. '''
        self.mergedict = mergedict
        self.states = states
        self.bins = {}

        self.bin_positions()

    def bin_positions(self):
        ''' Collapse equivalent states together. '''
        merged = [self.mergedict.get(state, state) for state in self.states]

        for index, state in enumerate(merged):
            if state not in self.bins:
                self.bins[state] = []
            self.bins[state].append(index)

    def shuffled(self, iterable=None, count=1, index=0):
        ''' Get a list of shuffled copies of an iterable

        iterable: Iterable to be shuffled. Falls back to list(range(len(states)))
        count: How many shuffles to generate.
        index: The lowest residue position if iterable is None. Primarily useful for explicit 0- or 1-indexing.
        '''

        iterable = list(range(index, index + len(self.states))) if iterable is None else iterable

        out = []

        bins = {}

        for _ in range(count):
            for state in self.bins: 
                bins[state] = self.bins[state][:]
                random.shuffle(bins[state])
            out.append([iterable[bins[self.mergedict.get(state, state)].pop()] for state in self.states])
        return out
        


def get_mergedict(inout=False, looptail=False):
    ''' Parse user-defined mergedict options into a mergedict 

    inout: Merge inside states [Ii] with outside states [Oo]
    looptail: Merge loop states [IO] with tail states [io] '''
    #TODO: optimize away all these annoyingly hardcoded lines
    #split both by tail/loop and in/out
    if inout and looptail: mergedict = {}
    #split by in/out but not tail/loop
    elif inout: mergedict = {'i':'I', 'o':'O'}
    #split by tail/loop but not in/out
    elif looptail: mergedict = {'O':'I', 'o':'i'}
    #everything outside a TMS is a loop
    else: mergedict = {'i':'I', 'O':'I', 'o':'I'}
    return mergedict

def old_main(args):
    ''' Deprecated legacy function '''
    mergedict = get_mergedict(args.split_surface, args.split_tails)

    f = open(args.o, 'w')

    fasta = Bio.SeqIO.parse(args.infile, 'fasta')

    bufsize = 100

    for seq in fasta:
        remaining = args.n
        if args.echo_first: f.write('>{}\n{}\n'.format(seq.name, seq.seq))

        for i in range(0, args.n, bufsize):
            if remaining < bufsize: count = remaining % bufsize
            else: count = bufsize

            sequences = old_shuffle_seq(seq, count=count, prefix=args.prefix, mergedict=mergedict)
            f.write('\n'.join(sequences))
            remaining -= count

class MSArray(object):
    ''' Helper class for fast columnwise rearrangement of MSAs '''
    def __init__(self, msa, arr=None):
        ''' Constructor for MSArray

        msa: Bio.Align.MultipleSeqAlignment object.
        arr: Char array representing MSA residues. Automatically computed if not given. '''
        self.msa = msa
        if arr is None: self.compute()
        else: self.arr = arr

    def compute(self):
        ''' Compute (Numpy) character array of MSA. '''
        self.arr = np.zeros((len(self.msa), self.msa.get_alignment_length()), dtype='<U1')
        for row, seqrecord in enumerate(self.msa):
            self.arr[row,:] = list(seqrecord.seq)

    def __getitem__(self, index):
        ''' Retrieve a slice from the array '''
        return self.arr.__getitem__(index)

    def reordered(self, order, prefix=None):
        ''' Get a new MSArray with a reodered MSA 

        order: 0-indexed new order of residues relative to the old one. 
        prefix: Prefix to add to sequence names. Liable to disappear in CLUSTAL format. '''
        newarr = np.zeros_like(self.arr)
        newarr = self.arr[:,list(order)]
        seqlist = []
        h = abs(hash(tuple(order)))
        for seqrecord, row in zip(self.msa, newarr):
            if prefix: seqid = '{}_{}_{:016x}'.format(seqrecord.id, prefix, h)
            else: seqid = seqrecord.id
            seqlist.append(Bio.SeqIO.SeqRecord(id=seqid, description='', seq=Bio.Seq.Seq(''.join(row))))
        newmsa = Bio.Align.MultipleSeqAlignment(seqlist)

        return MSArray(msa=newmsa, arr=newarr)

def main(args):
    validate_args(args)

    mergedict = get_mergedict(args.split_surface, args.split_tails)
    firstline = True
    bufsize = args.n

    if args.outfmt.startswith('numeric'):
        outmode = 'file'
        outfh = sys.stdout if args.outfile is None else open(args.outfile, 'w')

    elif args.outfile is None and args.outdir is None:
        outmode = 'file'
        outfh = sys.stdout

    elif args.outdir is None:
        outmode = 'file'
        outfh = open(args.outfile, 'w')

    elif args.outfile is None:
        outmode = 'dir'
        outdir = args.outdir
        outdir.mkdir(exist_ok=True)

    else: raise TypeError('At most one of args.outfile and args.outdir can be defined')

    #if numeric, begin preparing the output as a JSON object
    if args.outfmt.startswith('numeric'): outfh.write('[')

    if args.infmt == 'clustal': 
        for msa in Bio.AlignIO.parse(args.infile, 'clustal'):
            #no custom TMSs: slow consensus generation
            if args.indices is None: 
                states = []
                statecounts = [dict() for _ in range(msa.get_alignment_length())]
                firststates = None
                for seqrecord in msa:

                    if args.states is None: template, seqstates = get_pieces(seqrecord, mergedict=mergedict, topo=args.topo, indices=args.indices)
                    else: template, seqstates = list(seqrecord.seq), list(args.states)
                    if args.dump_states: 
                        outfh.write(format(Bio.SeqIO.SeqRecord(seq=Bio.Seq.Seq(''.join(seqstates)), id=str(pathlib.Path(args.infile).stem + '_states'), name='', description=''), 'fasta') + '\n')
                        break

                    for resi, state in enumerate(seqstates):
                        if state in statecounts[resi]: statecounts[resi][state] += 1
                        else: statecounts[resi][state] = 1
                    if firststates is None: firststates = seqstates

                if args.dump_states: continue

                for resi, counts in enumerate(statecounts):
                    #wherein there is only one option
                    if len(counts) == 1: states.append(list(counts)[0])
                    else:
                        ismax = [v == max(counts.values()) for v in counts.values()]
                        #wherein there is one obvious option
                        if sum(ismax) == 1:
                            [states.append(k) for k in counts if counts[k] == max(counts.values())]

                        else: #fall back on the first sequence I guess >:(
                            states.append(firststates[resi])
            #explicit custom TMSs: fast!
            else: 
                if args.states is None: template, states = get_pieces(msa[0], mergedict=mergedict, topo=args.topo, indices=parse_ranges(args.indices, offset=-1))
                else: template, states = list(msa[0].seq), list(args.states)
                if args.dump_states: 
                    outfh.write(format(Bio.SeqIO.SeqRecord(seq=Bio.Seq.Seq(''.join(states)), id=record.id + '_states', name='', description=''), 'fasta') + '\n')
                    continue

            #prepare the shuffler
            shuffler = StateShuffler(mergedict, states)
            msarr = MSArray(msa)
            remaining = args.n

            #shuffle and write as directed
            for i in range(0, args.n, bufsize):
                if remaining < bufsize: count = remaining % bufsize
                else: count = bufsize

                if args.outfmt.startswith('numeric'):
                    index = 1 if args.outfmt == 'numeric1' else 0
                    for numlist in shuffler.shuffled(count=count, index=index):
                        if firstline:
                            firstline = False
                            outfh.write('\n' + str(numlist))
                        else:
                            outfh.write(',\n' + str(numlist))
                else:
                    for order in shuffler.shuffled(count=count):
                        newmsarr = msarr.reordered(order, prefix=args.prefix)
                        if outmode == 'file': outfh.write(format(newmsarr.msa, args.outfmt))
                        else: 
                            h = hash(tuple(order))
                            ext = ''
                            if args.outfmt == 'clustal': ext = '.aln'
                            elif args.outfmt == 'fasta': ext = '.faa'
                            Bio.AlignIO.write(newmsarr.msa, outdir.joinpath('{}_{:016X}{}'.format(args.prefix, abs(h), ext)), args.outfmt)
                remaining -= args.n

    elif args.infmt == 'fasta':
        for record in Bio.SeqIO.parse(args.infile, 'fasta'):
            remaining = args.n

            #add the original if requested (unless numeric format was requested)
            if args.echo_first and args.outfmt != 'numeric': fh.write(format(record, args.outfmt) + '\n')

            #prepare the shuffler
            if args.indices is None:
                topo, indices = hmmtop(record)
            else: 
                topo = args.topo
                indices = parse_ranges(args.indices, offset=-1)

            if args.states is None: template, states = get_pieces(record, mergedict=mergedict, topo=topo, indices=indices)
            else: template, states = list(record.seq), list(args.states)
            if args.dump_states: 
                outfh.write(format(Bio.SeqIO.SeqRecord(seq=Bio.Seq.Seq(''.join(states)), id=record.id + '_states', name='', description=''), 'fasta') + '\n')
                continue
            shuffler = StateShuffler(mergedict, states)

            #shuffle and write as directed
            for i in range(0, args.n, bufsize):
                if remaining < bufsize: count = remaining % bufsize
                else: count = bufsize

                if args.outfmt.startswith('numeric'):
                    index = 1 if args.outfmt == 'numeric1' else 0
                    for numlist in shuffler.shuffled(count=count, index=index):
                        if firstline:
                            firstline = False
                            outfh.write('\n' + str(numlist))
                        else:
                            outfh.write(',\n' + str(numlist))
                else:
                    for reslist in shuffler.shuffled(record.seq, count=count): 
                        seq = ''.join(reslist)
                        newrecord = Bio.SeqIO.SeqRecord(id='{}_{:016X}'.format(args.prefix, abs(hash(seq))), description='', seq=Bio.Seq.Seq(seq))
                        if outmode == 'file':
                            outfh.write(format(newrecord, args.outfmt))
                        else:
                            h = abs(hash(seq))
                            Bio.SeqIO.write(newrecord, outdir.joinpath('{}_{:016X}.faa'.format(args.prefix, abs(h))), args.outfmt)
                remaining -= args.n
            
    else: raise NotImplementedError

    #finish the JSON
    if args.outfmt.startswith('numeric'): outfh.write('\n]')

def validate_args(args):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--echo-first', action='store_true', help='Echo seed sequence in output as first sequence')
    parser.add_argument('-n', default=1, type=int, help='Number of shuffles (default:1)')
    outparser = parser.add_mutually_exclusive_group()
    outparser.add_argument('-o', '--outfile', type=pathlib.Path, default=None, help='Where to write the shuffled sequences (default:stdout)')
    outparser.add_argument('--outdir', type=pathlib.Path, default=None, help='Where to write the shuffled sequences. --outfmt numeric reverts to --outfile output.')
    parser.add_argument('-p', '--prefix', default='shuffled', help='Sequence label prefix (default:shuffled)')
    parser.add_argument('-s', '--split-surface', action='store_true', help='Split I/i from O/o when shuffling to keep loop topology')
    parser.add_argument('-t', '--split-tails', action='store_true', help='Split I from i and O from o when shuffling to retain some information on distance from the membrane')
    parser.add_argument('--infmt', default='fasta', help='Input format. If clustal, this script shuffles relevant positions once for each MSA.')
    parser.add_argument('--outfmt', default='fasta', help='Output format. Tested formats: fasta, numeric (0-indexed numeric), numeric1 (1-indexed numeric), clustal (if infmt is clustal). (default:fasta)')

    #customtms = parser.add_mutually_exclusive_group()
    parser_customtms = parser.add_argument_group('arguments for working with custom TMSs')
    parser_customtms.add_argument('--topo', type=int, default=None, help='Initial residue environment. Accepted values: -1: outer loop, 0: membrane, 1: cytosolic.')
    parser_customtms.add_argument('--indices', help='TMS coordinates as comma-separated dash-delimited ranges, e.g. 20-40,60-80,200-230')

    parser_customtms.add_argument('--states', help='Custom states string. Must match input sequence in length')

    parser.add_argument('--buffer', type=int, default=None, help='Size of write buffer in shuffles. (default:100)')

    parser.add_argument('--dump-states', action='store_true', help='Dump states instead of shuffling')

    parser.add_argument('infile', nargs='?', default='/dev/stdin', help='File to read in (default:stdin)')

    args = parser.parse_args()

    if sum([(args.topo is not None or args.indices is not None), args.states is not None]) >= 2:
        parser.error('--topo/--indices and --states are mutually exclusive.')
        exit(1)


    main(args)

