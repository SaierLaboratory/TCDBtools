import subprocess
import re
import io
import os

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm

from Bio import SeqIO, Seq, Align, AlignIO
from kevtools_common import badnan
from scipy.stats import gaussian_kde
import argparse
import sys
import tempfile

from . import indexing

def info(*things):
    print('[INFO]', *things, file=sys.stderr)
def warn(*things):
    print('[WARNING]', *things, file=sys.stderr)
def error(*things):
    print('[ERROR]', *things, file=sys.stderr)
    exit(1)

def apply_fragment(fragseq, fullseq, arr):

    aligner = Align.PairwiseAligner(mode='local', match_score=2, mismatch_score=-1, gap_score=0)
    alignments = aligner.align(str(fragseq.seq).replace('-', 'X'), str(fullseq.seq).replace('-', 'X'))

    first = alignments[0]
    fragresi = []
    fullresi = []
    for start, end in first.aligned[0]: fragresi.extend(range(start, end))
    for start, end in first.aligned[1]: fullresi.extend(range(start, end))

    shortfrag = str(fragseq.seq).replace('-', '')
    shortfull = str(fullseq.seq).replace('-', '')

    #STEP 1: SANITY CHECKS
    #if fragment isn't a subsequence of full, crash
    if len(shortfrag) > len(shortfull): raise IndexError('fragseq ({}) is longer than fullseq ({})'.format(fragseq.id, fullseq.id))
    elif (fragresi[0] > 0) or (fragresi[-1] < (len(fragseq) - 1)): raise IndexError('fragseq ({}) is not a subsequence of fullseq ({})'.format(fragseq.id, fullseq.id))
    elif len(fragresi) != len(fullresi): raise IndexError('aligned residues for fragseq ({}) don\'t match aligned residues for fullseq ({})'.format(fragseq.id, fullseq.id))
    #if arr is longer than fullseq, crash
    if len(arr) > len(fullseq): raise IndexError('fullseq ({}) is longer than data array'.format(fullseq.id))
    #if fullseq references resi outside arr range, crash
    if max(fullresi) >= len(arr) or min(fullresi) < 0: raise IndexError('fullseq ({}) contains residues not in data array'.format(fullseq.id))

    newarr = np.nan * np.zeros(len(fragseq), dtype=arr.dtype)
    #STEP 2: COPY STUFF OVER
    for fragi, fulli in zip(fragresi, fullresi):
        newarr[fragi] = arr[fulli]
    return newarr



    #newarr = np.zeros(len(fragseq), dtype=arr.dtype)
    #fragi = fragresi[0]
    #indexi = 0
    #for c in fragseq:
    #    if c in 'ACDEFGHIKLMNPQRSTVWY': 
    #        fragi[indexi] = arr[fulli[indexi]]
    #        indexi += 1

    fn1 = tempfile.NamedTemporaryFile('w')
    SeqIO.write([fragseq], fn1, 'fasta')
    fn1.flush()
    fn2 = tempfile.NamedTemporaryFile('w')
    SeqIO.write([fullseq], fn2, 'fasta')
    fn2.flush()

    cmd = ['blastp', '-query', fn1.name, '-subject', fn2.name, '-outfmt', '6 pident qstart qend qlen sstart send slen qseq sseq']
    out = subprocess.check_output(cmd).decode('utf-8')

    if not out.strip(): error('Fragment ({}) and full sequence ({}) are unalignable'.format(fragseq.id, fullseq.id))
    else: 
        #frag, full = out.strip().split('\n')[0].split('\t')
        firstmatch = out.strip().split('\n')[0].split('\t')
        if firstmatch[0] != '100.000': warn('Fragment ({}) and full sequence ({}) are not identical'.format(fragseq.id, fullseq.id))

        if (firstmatch[1] != '1') or (firstmatch[2] != firstmatch[3]): warn('Fragment ({}) is not fully covered by full sequence ({})'.format(fragseq.id, fullseq.id))

        newarr = np.zeros(0)

        ### insertions, deletions

        if ('-' in firstmatch[-2]) or ('-' in firstmatch[-1]): raise NotImplementedError('Please use exact sequences without indels')

        else:
            sstart = int(firstmatch[4])
            send = int(firstmatch[5])
            newarr = arr[sstart-1:send]

        #deletions
        #if '-' in firstmatch[-2]: 
        #    gaps = set()
        #    for i, c in firstmatch[-2]:
        #        if c == '-': gaps.add(i)
        #info(newarr)
        return newarr

def close_gaps(shiftedindices):
    workingcopy = shiftedindices.copy()
    nans = np.nonzero(np.isnan(shiftedindices))[0]

    for startindex in nans:
        if startindex > 0 and shiftedindices[startindex - 1] == 0:
            workingcopy[startindex] = 0
            return close_gaps(workingcopy)
        elif startindex < (len(shiftedindices) - 1) and shiftedindices[startindex + 1] == 0:
            workingcopy[startindex] = 0
            return close_gaps(workingcopy)

    workingcopy[np.isnan(workingcopy)] = 1
    return workingcopy

def shannon(seq, window=20):
    #FIXME: Implement the O(n) algo instead of this O(n**2)-like one

    shortseq, nans = badnan.collapse(seq, '-')
    vals = np.zeros(len(shortseq)-window+1)
    for i in range(len(shortseq)-window+1):
        compdict = {}
        for j in range(window):
            resn = shortseq[i+j]
            if resn in compdict: compdict[resn] += 1
            else: compdict[resn] = 1

        total = sum(compdict.values())
        for resn in compdict:
            p = compdict[resn]/total
            vals[i] += - p * np.log2(p)

    entropies = vals
    entropies = np.hstack([np.zeros(window//2), entropies, np.zeros(window - window//2)])
    entropies = badnan.uncollapse(entropies, nans, '-')

    entropies[:window//2] = np.nan
    entropies[-(window-window//2):] = np.nan

    entropies = entropies[window//2:-(window-window//2)]
    return entropies

def str2seq(text):
    cleantext = re.sub('[^-A-Z]', '', text)
    seqid = 'seq{:016x}'.format(abs(hash(cleantext)))
    seq = Seq.Seq(cleantext)
    return SeqIO.SeqRecord(seq, id=seqid, name=seqid, description=seqid)


def update_lims(oldxlim, oldylim, newlims=None):
    if newlims is None: return oldxlim, oldylim

    if oldxlim is None: xmin = xmax = None
    else: xmin, xmax = oldxlim

    if oldylim is None: ymin = ymax = None
    else: ymin, ymax = oldylim

    if newlims[0] is None: pass
    else:
        if newlims[0][0] is not None: 
            if xmin is None or (newlims[0][0] < xmin): xmin = newlims[0][0]
        if newlims[0][1] is not None: 
            if xmax is None or (newlims[0][1] > xmax): xmax = newlims[0][1]
    if newlims[1] is None: pass
    else:
        if newlims[1][0] is not None: 
            if ymin is None or (newlims[1][0] < ymin): ymin = newlims[1][0]
        if newlims[1][1] is not None: 
            if ymax is None or (newlims[1][1] > ymax): ymax = newlims[1][1]

    return (xmin, xmax), (ymin, ymax)

class Kernels(object):
    @staticmethod
    def moving_average(window):
        return (1/window) * np.ones(window)

    @staticmethod
    def cosine(window):
        k = (1 - np.cos(2*np.pi/window*np.arange(window)))
        return k / np.sum(k)

    @staticmethod
    def triangular(window):
        kernel = (1 - 2/window * np.abs(np.arange(window)-window//2)) / window * 2
        return kernel

    @staticmethod
    def tukey(window, a=0.5):
        kernel = np.zeros(window)
        for i in np.arange(-window//2, window-window//2):
            if 0 <= i < a*window/2: kernel[i] = 0.5 * (1 - np.cos(2*np.pi*i/a/window))
            elif a*window/2 <= i <= window/2: kernel[i] = 1
            else: kernel[i] = 0.5 * (1 - np.cos(2*np.pi*(window-i)/a/window))
        return kernel/sum(kernel)

    @staticmethod
    def get_kernel(window, name=None):
        if name is None: return Kernels.moving_average(window)
        elif name.startswith('moving') or name.startswith('flat') or name.startswith('rect') or name in ('-', '_'): return Kernels.moving_average(window)
        elif name.startswith('hann') or name.startswith('cos'): return Kernels.cosine(window)
        elif name.startswith('triang'): return Kernels.triangular(window)
        elif name.startswith('tukey'): return Kernels.tukey(window)
        else: error('Kernel not found: "{}"'.format(name))

class Tables(object):
    hydropathy = {'G':-0.400, 
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
     'B':-3.500, 
     'J':4.150, 
     'Z':-3.500, 
    }
    kyte_dolittle_hydrophobicity = {'G':-0.400, 
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
     'B':-3.500, 
     'J':4.150, 
     'Z':-3.500, 
    }
    eisenberg_hydrophobicity = {'G':0.48, 
     'I':1.38, 
     'S':-0.18, 
     'Q':-0.85, 
     'E':-0.74, 
     'A':0.62, 
     'M':0.64, 
     'T':-0.05, 
     'Y':0.26, 
     'H':-0.4, 
     'V':1.08, 
     'F':1.19, 
     'C':0.29, 
     'W':0.81, 
     'K':-1.5, 
     'L':1.06, 
     'P':0.12, 
     'N':-0.78, 
     'D':-0.90, 
     'R':-2.53, 
     'B':-0.84, 
     'J':1.22, 
     'Z':-0.8,
    }

    conformational = {
        'A':2,
        'B':4,
        'C':4,
        'D':4,
        'E':5,
        'F':4,
        'G':3,
        'H':4,
        'I':4,
        'J':4,
        'K':6,
        'L':4,
        'M':5,
        'N':4,
        'O':6,
        'P':1,
        'Q':5,
        'R':6,
        'S':4,
        'T':4,
        'U':4,
        'V':3,
        'W':4,
        'Y':5,
        'Z':5,
    }

    charge = {
        'A':0.0,
        'B':-0.5,
        'C':+0.1,
        'D':-1.0,
        'E':-1.0,
        'F':0.0,
        'G':0.0,
        'H':+0.1,
        'I':0.0,
        'J':-0.5,
        'K':+1.0,
        'L':0.0,
        'M':0.0,
        'N':0.0,
        'O':0.0,
        'P':0.0,
        'Q':0.0,
        'R':+1.0,
        'S':0.0,
        'T':0.0,
        'U':0.0,
        'V':0.0,
        'W':0.0,
        'Y':+0.001,
        'Z':0.0,
    }

    cfhelix = {
        "A":1.42,
        "C":0.70,
        "D":1.01,
        "E":1.51,
        "F":1.13,
        "G":0.57,
        "H":1.00,
        "I":1.08,
        "K":1.16,
        "L":1.21,
        "M":1.45,
        "N":0.67,
        "P":0.57,
        "Q":1.11,
        "R":0.98,
        "S":0.77,
        "T":0.83,
        "V":1.06,
        "W":1.08,
        "Y":0.69,
    }
    cfsheet = {
        "A":0.83,
        "C":1.19,
        "D":0.54,
        "E":0.37,
        "F":1.38,
        "G":0.75,
        "H":0.87,
        "I":1.60,
        "K":0.74,
        "L":1.30,
        "M":1.05,
        "N":0.89,
        "P":0.55,
        "Q":1.10,
        "R":0.93,
        "S":0.75,
        "T":1.19,
        "V":1.70,
        "W":1.37,
        "Y":1.47,
    }



    @staticmethod
    def get_table(name=None):
        if name is None: return Tables.kyte_dolittle_hydrophobicity
        elif name.startswith('hydro'): return Tables.kyte_dolittle_hydrophobicity
        elif name.startswith('amphi'): return Tables.eisenberg_hydrophobicity
        elif name.startswith('conform'): return Tables.conformational
        elif name.startswith('charge'): return Tables.charge
        else: raise NotImplementedError('Table not found: "{}"'.format(name))

def lookup(sequence, table):

    return np.array([table.get(r, np.nan) for r in sequence])

def nanconvolve(arr1, arr2, mode='same', pad=True, ignore_oob=True): 
    if len(arr1) < len(arr2): 
        longer = arr2
        shorter = arr1
    else: 
        longer = arr1
        shorter = arr2

    nanless, positions = badnan.collapse(longer)

    #maxlen = max(len(arr1), len(arr2))
    minlen = min(len(arr1), len(arr2))

    result = badnan.uncollapse(np.convolve(nanless, shorter, mode=mode), positions, ignore_oob=ignore_oob)

    #if pad: result = np.pad(result, (minlen//2, minlen - minlen//2 - 1), mode='constant', constant_values=np.nan)
    if pad:
        result[:minlen//2] = np.nan
        result[-minlen//2:] = np.nan
    else:
        result = result[minlen//2:-minlen//2]

    return result

class BaseEntity(object):
    label = None
    color = None
    facecolor = None
    edgecolor = None
    textcolor = None
    linewidth = None
    
    def __init__(self, **kwargs):
        pass

    @staticmethod
    def compute(seq, **kwargs): raise NotImplementedError

    @staticmethod
    def compute_msa(msa, **kwargs): raise NotImplementedError

    def plot(self, ax=None): raise NotImplementedError

    def get_bounding_box(self): raise NotImplementedError

    def set_color(self, color):
        self.color = color
        self.facecolor = color
        self.edgecolor = color

    def set_edgecolor(self, color):
        self.edgecolor = color

    def set_facecolor(self, color):
        self.facecolor = color

    def apply_styleargs(self, **kwargs):
        if "style" in kwargs: self.color = kwargs["style"]

        if "color" in kwargs: self.color = kwargs["color"]
        if "alpha" in kwargs: self.alpha = kwargs["alpha"]

        if "ec" in kwargs: self.edgecolor = kwargs["ec"]
        if "edgecolor" in kwargs: self.edgecolor = kwargs["edgecolor"]
        if "fc" in kwargs: self.facecolor = kwargs["fc"]
        if "facecolor" in kwargs: self.facecolor = kwargs["facecolor"]

        if type(self.facecolor) is str and kwargs.get("alpha") is None:
            if (self.facecolor.startswith("#") and len(self.facecolor) == 9) or (self.facecolor.startswith("0x") and len(self.facecolor) == 10):
                self.alpha = None

        if "lw" in kwargs: self.linewidth = kwargs["lw"]
        if "linewidth" in kwargs: self.linewidth = kwargs["linewidth"]

        if "textcolor" in kwargs: self.textcolor = kwargs.get('textcolor', None)


        if self.color is not None:
            if self.edgecolor is None: self.edgecolor = self.color
            if self.facecolor is None: self.facecolor = self.color
        if type(self.edgecolor) is int: self.edgecolor = get_darkcolor(self.edgecolor)
        if type(self.facecolor) is int: self.facecolor = get_lightcolor(self.facecolor)


class BaseCurve(BaseEntity):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.X = kwargs.get('X', np.zeros(0))
        self.Y = kwargs.get('Y', np.zeros(0))

        self.label = kwargs.get('label', 'BaseCurve')

        self.apply_styleargs(**kwargs)

        self.record = None


    def plot(self, ax=None): ax.plot(self.X, self.Y, label=self.label, lw=self.linewidth, color=self.edgecolor)

class Hydropathy(BaseCurve):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.X = np.zeros(0)
        self.Y = np.zeros(0)
        self.window = kwargs.get('window', 19)
        self.kernel = kwargs.get('kernel', None)
        self.table = kwargs.get('table', None)

        self.label = 'Hydropathy'

        self.apply_styleargs(**kwargs)


    @staticmethod
    def compute(seq, **kwargs):

        window = kwargs.get('window', 19)
        kernel = kwargs.get('kernel', None)
        table = kwargs.get('table', None)
        fragment = kwargs.get('fragment', None)

        X = np.arange(len(seq))
        rawY = lookup(seq.seq, Tables.get_table(table))
        Y = nanconvolve(rawY, Kernels.get_kernel(window=window, name=kernel))
        hydro = Hydropathy()
        hydro.X = X
        hydro.Y = Y
        hydro.record = seq

        if fragment: 
            hydro.Y = apply_fragment(fragment, seq, hydro.Y)
            hydro.X = np.arange(len(hydro.Y))

        hydro.apply_styleargs(**kwargs)

        return hydro

    def compute_msa(msa, **kwargs):
        window = kwargs.get('window', 19)
        kernel = kwargs.get('kernel', None)
        table = Tables.get_table(kwargs.get('table', None))
        #fragment = kwargs.get('fragment', None) #not implemented yet

        positions = np.zeros(msa.get_alignment_length())
        occupancies = np.zeros(msa.get_alignment_length())

        for seq in msa:
            for i, resn in enumerate(seq):
                if resn in table:
                    if table[resn]:
                        positions[i] += table[resn]
                        occupancies[i] += 1

        positions /= occupancies
        X = np.arange(msa.get_alignment_length())
        rawY = positions / occupancies
        collY, gaps = indexing.collapse(positions)
        convY = nanconvolve(collY, Kernels.get_kernel(name=kernel, window=window))
        Y = indexing.uncollapse(list(convY), gaps)
        if window > 1:
            Y[:window//2] = np.nan
            Y[-window//2:] = np.nan
        
        hydro = Hydropathy()
        hydro.apply_styleargs(**kwargs)

        hydro.X = X
        hydro.Y = Y

        #do fragment processing?

        return hydro

    def plot(self, ax=None):
        if ax is None:
            ax = plt.figure().gca()
            ax.axhline(0, color='black', lw=1)
            ax.set_ylim([-3, 3])

        super().plot(ax)
        #ax.plot(self.X, self.Y, color=self.edgecolor, label=self.label, linewidth=self.linewidth)

        if ax is None:
            plt.show()

    def get_bounding_box(self):
        xlim = np.min(self.X), np.max(self.X) + 1
        ylim = -3, 3
        return xlim, ylim

class Charge(Hydropathy):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.table = kwargs.get('table', 'charge')
        self.label = 'Charge'
    @staticmethod
    def compute(seq, **kwargs):
        window = kwargs.get('window', 21)
        kernel = kwargs.get('kernel', None)
        table = kwargs.get('table', 'charge')
        fragment = kwargs.get('fragment', None)
        

        X = np.arange(len(seq))
        rawY = lookup(seq.seq, Tables.get_table(table))
        Y = nanconvolve(rawY, Kernels.get_kernel(window=window, name=kernel))
        Y[:window//2] = np.nan
        Y[-window//2:] = np.nan
        charge = Charge()
        charge.apply_styleargs(**kwargs)
        charge.X = X
        charge.Y = Y
        charge.record = seq

        if fragment: 
            charge.Y = apply_fragment(fragment, seq, charge.Y)
            charge.X = np.arange(len(charge.Y))

        return charge

class Entropy(Hydropathy):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.label = 'Compositional entropy'

    @staticmethod
    def compute(seq, **kwargs):
        window = kwargs.get('window', 21)
        kernel = kwargs.get('kernel', None)
        fragment = kwargs.get('fragment', None)

        rawY = shannon(seq.seq, window)
        
        if kernel == 'flat':
            Y = rawY
            X = np.arange(len(Y)) + window//2
        else:
            #FIXME: double smoothing when two degrees of freedom vanish
            Y = nanconvolve(rawY, Kernels.get_kernel(window=window, name=kernel))
            Y[:window//2] = np.nan
            Y[-window//2:] = np.nan
            X = np.arange(len(Y)) + window//2


        entropy = Entropy()
        entropy.apply_styleargs(**kwargs)
        entropy.X = X
        entropy.Y = Y
        entropy.record = seq

        if fragment: 
            entropy.Y = apply_fragment(fragment, seq, entropy.Y)
            entropy.X = np.arange(len(entropy.Y))

        return entropy

    def get_bounding_box(self):
        xlim = np.min(self.X), np.max(self.X)
        ylim = 2, 5
        return xlim, ylim
        
class Conformational(Hydropathy):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.X = np.zeros(0)
        self.Y = np.zeros(0)
        self.window = kwargs.get('window', 19)
        self.kernel = kwargs.get('kernel', None)
        self.table = kwargs.get('table', 'conformational')

        self.label = 'Conformational entropy'

        self.linewidth = None

    @staticmethod
    def compute(seq, **kwargs):
        window = kwargs.get('window', 19)
        kernel = kwargs.get('kernel', None)
        table = kwargs.get('table', 'conformational')
        fragment = kwargs.get('fragment', None)

        X = np.arange(len(seq))
        rawY = lookup(seq.seq, Tables.get_table(table))
        Y = nanconvolve(rawY, Kernels.get_kernel(window=window, name=kernel))
        hydro = Conformational()
        hydro.apply_styleargs(**kwargs)
        hydro.X = X
        hydro.Y = Y
        hydro.record = seq

        if fragment: 
            hydro.Y = apply_fragment(fragment, seq, hydro.Y)
            hydro.X = np.arange(len(hydro.Y))

        return hydro

    def get_bounding_box(self):
        xlim = np.min(self.X), np.max(self.X)
        ylim = 2, 5
        return xlim, ylim

class Identity(Hydropathy):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.X = np.zeros(0)
        self.Y = np.zeros(0)
        self.window = kwargs.get('window', 19)
        self.kernel = kwargs.get('kernel', 19)

        self.label = 'Identity'


    @staticmethod
    def compute(seq1, seq2, **kwargs):
        window = kwargs.get('window', 19)
        kernel = kwargs.get('kernel', None)
        fragment = kwargs.get('fragment', None)

        X = np.arange(min(len(seq1), len(seq2)))
        #FIXME: handle double gaps and such
        rawY = np.array([1 if (r1 == r2) else 0 for r1, r2 in zip(seq1.seq, seq2.seq)])
        Y = nanconvolve(rawY, Kernels.get_kernel(window=window, name=kernel))
        hydro = Identity()
        hydro.apply_styleargs(**kwargs)
        hydro.X = X
        hydro.Y = Y
        hydro.record = seq1
        return hydro

    def get_bounding_box(self):
        xlim = np.min(self.X), np.max(self.X)
        ylim = 0, 1
        return xlim, ylim

class Occupancy(Hydropathy):
    ''' For MSAs and profiles '''

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.X = np.zeros(0)
        self.Y = np.zeros(0)
        self.window = kwargs.get('window', 1)
        self.kernel = kwargs.get('kernel')

        self.ymin = kwargs.get('ymin', -3)
        self.ymax = kwargs.get('ymax', -2)

        self.edgecolor = 'y'


        self.label = 'Occupancy'

    @staticmethod
    def compute(seq, **kwargs):
        window = kwargs.get('window', 1)
        kernel = kwargs.get('kernel')

        ymin = kwargs.get('ymin', -3)
        ymax = kwargs.get('ymax', -2)

        X = np.arange(len(seq.seq))
        rawY = np.ones(len(seq.seq))

        for i, resn in enumerate(seq): 
            if resn == '-': rawY[i] = 0

        if window == 1: Y = rawY
        else:
            Y = nanconvolve(rawY, Kernels.get_kernel(window, kernel))

        Y = (ymax - ymin) * Y + ymin

        obj = Occupancy()
        obj.apply_styleargs(**kwargs)
        obj.X = X
        obj.Y = Y
        obj.record = seq
        return obj

    @staticmethod
    def compute_msa(msa, **kwargs):
        window = kwargs.get('window', 1)
        kernel = kwargs.get('kernel')

        ymin = kwargs.get('ymin', -3)
        ymax = kwargs.get('ymax', -2)

        X = np.arange(msa.get_alignment_length())
        Y = np.zeros(msa.get_alignment_length())
        for record in msa:
            for alni, resn in enumerate(record.seq):
                if resn in set(('-', '_', 'X')): pass
                else: Y[alni] += 1

        Y /= len(msa)

        if window == 1: pass
        else:
            Y = nanconvolve(Y, Kernels.get_kernel(window, kernel))

        Y = (ymax - ymin) * Y + ymin

        obj = Occupancy()
        obj.apply_styleargs(**kwargs)
        obj.X = X
        obj.Y = Y
        return obj

    def get_bounding_box(self):
        xlim = np.min(self.X), np.max(self.X)
        ylim = 0, 1
        return xlim, ylim

class AveHASimilarity(Hydropathy):
    ''' AveHAS-style similarity, i.e. -stdev(*hydropathy) '''

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.ymin = kwargs.get('ymin', -3)
        self.ymax = kwargs.get('ymax', -1.5)

        self.gap_wt = kwargs.get('gap_wt', 1.)
        self.hydro_wt = kwargs.get('hydro_wt', 1.)

        self.label = 'Similarity'
    @staticmethod
    def compute(seq, **kwargs):
        #makes no sense to do this calculation, but oh well

        obj = AveHASimilarity()
        ymax = kwargs.get('ymax', obj.ymax)

        X = np.arange(len(seq))
        Y = np.ones(len(seq)) * ymax

        obj.apply_styleargs(**kwargs)
        obj.X = X
        obj.Y = Y
        obj.record = seq
        return obj

    @staticmethod
    def compute_msa(msa, **kwargs):

        table = kwargs.get('table', 'hydro')
        #need more robust type checking for this sort of thing, probably
        if not isinstance(table, dict): table = Tables.get_table(table)

        kernel = kwargs.get('kernel')
        window = kwargs.get('window', 11)


        obj = AveHASimilarity(window=window, kernel=kernel, table=table)

        gap_wt = kwargs.get('gap_wt', 1.)
        hydro_wt = kwargs.get('hydro_wt', 1.)

        ymin = kwargs.get('ymin', obj.ymin)
        ymax = kwargs.get('ymax', obj.ymax)

        gap_counts = np.zeros(msa.get_alignment_length())
        gap_costs = np.zeros(msa.get_alignment_length())
        hydro_stdevs = np.zeros(msa.get_alignment_length())

        for i in range(msa.get_alignment_length()):
            position = []
            gaps = 0
            for seq in msa:
                resn = seq[i]
                if resn == '-': gaps += 1
                else: position.append(table.get(resn, np.nan))
            gap_costs[i] = (len(msa) - gaps - 1) / (len(msa) - 1)

            score = -np.nanstd(position)

            hydro_stdevs[i] = score
        minstd, maxstd = np.min(hydro_stdevs), np.max(hydro_stdevs)

        norm_hydro_stdevs = (hydro_stdevs - minstd) / (maxstd - minstd)
        scores = norm_hydro_stdevs**hydro_wt * gap_costs**gap_wt

        collscores, gaps = indexing.collapse(scores)
        convscores = nanconvolve(collscores, Kernels.get_kernel(name=kernel, window=window))
        Y = indexing.uncollapse(convscores, gaps)
        if window > 1:
            Y[:window//2] = np.nan
            Y[-window//2:] = np.nan

        Y = Y * (ymax - ymin) + ymin
        X = np.arange(len(Y))

        obj.apply_styleargs(**kwargs)
        obj.X = X
        obj.Y = Y

        return obj


class Amphipathicity(Hydropathy):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.label = 'Amphipathicity'

    @staticmethod
    def compute(seq, **kwargs):
        kernelname = kernel = kwargs.get('kernel') #multiplied with rotation kernel as windowing 
        window = kwargs.get('window', 19)
        table = kwargs.get('table', 'amphi')
        fragment = kwargs.get('fragment')

        angle = kwargs.get('angle', 100)
        offsets = np.arange(window) - window//2
        rawwindow = np.cos(360/angle * 2 * np.pi * offsets)

        scale = kwargs.get('scale', 2.5) #I guess it looks better with this?

        #if kernel is None or kernel.startswith('flat'): kernel = rawwindow / window
        #else:
        #    kernel = Kernels.get_kernel(window=window, name=kernel) * rawwindow / window * Kernels.get_kernel(window=window, name='gaussian')
        kernel = Kernels.get_kernel(window=window, name=kernel)
        kernel /= sum(kernel)

        rawY = lookup(seq.seq, Tables.get_table(table))
        collY, gaps = indexing.collapse(rawY)
        cosines = collY * np.cos(np.pi / 180 * angle * np.arange(len(collY)))
        sines = collY * np.sin(np.pi / 180 * angle * np.arange(len(collY)))
        convY = np.sqrt(nanconvolve(sines, kernel)**2 + nanconvolve(cosines, kernel)**2) * scale
        Y = indexing.uncollapse(convY, gaps)
        if window > 1:
            Y[:window//2] = np.nan
            Y[-window//2:] = np.nan

        #X = np.arange(len(seq))
        #FIXME: Failure to expand properly with indexing.uncollapse
        X = np.arange(len(Y))

        obj = Amphipathicity()
        obj.apply_styleargs(**kwargs)
        obj.X = X
        obj.Y = np.abs(Y)
        obj.record = seq

        if fragment:
            obj.Y = apply_fragment(fragment, seq, obj.Y)
            obj.X = np.arange(len(obj.Y))

        return obj

    @staticmethod
    def compute_msa(msa, **kwargs):
        kernelname = kernel = kwargs.get('kernel') #multiplied with rotation kernel as windowing 
        window = kwargs.get('window', 19)
        table = kwargs.get('table', 'amphi')
        fragment = kwargs.get('fragment')

        angle = kwargs.get('angle', 100)
        offsets = np.arange(window) - window//2
        rawwindow = np.cos(360/angle * 2 * np.pi * offsets)

        scale = kwargs.get('scale', 2.5) #I guess it looks better with this?

        kernel = Kernels.get_kernel(window=window, name=kernel)
        kernel /= sum(kernel)

        rawY = np.nanmean(np.stack([lookup(seq.seq, Tables.get_table(table)) for seq in msa]), axis=0)
        length = msa.get_alignment_length()
        cosines = rawY * np.cos(np.arange(length) * np.pi / 180 * angle)
        sines = rawY * np.sin(np.arange(length) * np.pi / 180 * angle)
        occupancies = np.mean(np.stack([[resn != "-" for resn in seq.seq] for seq in msa]), axis=0)
        
        avgcosines = nanconvolve(occupancies * cosines, kernel)
        avgsines = nanconvolve(occupancies * sines, kernel)
        amphi = (avgcosines**2 + sines**2) ** .5
        amphi *= 5/2
        amphi /= 4.18

        obj = Amphipathicity() #FIXME: Compute this
        obj.apply_styleargs(**kwargs)
        obj.X = np.arange(len(amphi))
        obj.Y = amphi
        #obj.msa = msa
        return obj


class HMMTOP(BaseEntity):
    def __init__(self, spans=None, **kwargs):

        self.spans = [] if spans is None else spans
        self.centers = kwargs.get('centers')
        self.length = 0
        self.height = 1

        self.facecolor = "orange"
        self.alpha = 0.25
        self.apply_styleargs(**kwargs)

        self.mode = kwargs.get('mode')

        #fraction of peak required to assign as TMS for density MSA mode (default: 50%)
        self.threshold = kwargs.get('threshold', 0.5)

        #maybe transfer to usual HMMTOP plots?
        self.ymin = kwargs.get('ymin')
        self.ymax = kwargs.get('ymax')

    @staticmethod
    def parse(fh):
        spans = []
        for l in fh:
            sl = re.findall('(?:(?: *[0-9]+)+)$', l.strip())
            if not sl: continue
            else: sl = [int(x) for x in sl[0].strip().split()]
            if len(sl) % 2: sl = sl[1:]
            spans.extend([sl[i:i+2] for i in range(0, len(sl), 2)])
        return spans

    @staticmethod
    def parse_hmmtop(topout):
        if not topout: topout = '>HP: 0 sequence OUT   0       \n'
        return [[int(x) for x in re.findall('((?:[0-9]+\s*)+)$', topline)[0].split()] for topline in topout.strip().split('\n')]

    @staticmethod
    def compute(seq, **kwargs):
        if kwargs.get('notms'): return HMMTOP([])

        fragment = kwargs.get('fragment')

        p = subprocess.Popen(['hmmtop', '-if=--', '-of=--', '-pi=spred', '-sf=FAS', '-is=pseudo'], 
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        out, err = p.communicate(input='>sequence\n{}'.format(seq.seq).encode('utf-8'))
        print(err.decode('utf-8').strip(), file=sys.stderr)
        out = out.decode('utf-8')

        #indices = [int(x) for x in re.findall('((?:[0-9]+\s*)+)$', out.decode('utf-8'))[0].split()]
        indices = HMMTOP.parse_hmmtop(out)[0]
        offset = 0
        partial2full = {}
        for resi, resn in enumerate(seq.seq):
            if resn not in 'ACDEFGHIKLMNPQRSTVWY': offset += 1
            partial2full[resi+1] = resi + offset

        if len(indices) > 1:
            spans = [[partial2full[x] for x in indices[i:i+2]] for i in range(1, len(indices), 2)]
        else: spans = []

        if fragment:
            arr = np.zeros(len(seq), np.bool8)
            unrolled = indexing.unroll_ranges(spans)
            arr[[x-1 for x in unrolled]] = 1
            #print(indexing.roll_ranges(np.nonzero(arr)[0]))
            shiftedindices = close_gaps(apply_fragment(fragment, seq, arr))
            #print(indexing.roll_ranges(np.nonzero(shiftedindices == 1)[0]))
            spans = indexing.roll_ranges(np.nonzero(shiftedindices == 1)[0])
            #print(indexing.roll_ranges(np.nonzero(shiftedindices)))
            #newspans = indexing.roll_ranges(shiftedindices)
            #print(newspans)

        obj = HMMTOP(spans)
        obj.length = len(seq)
        obj.apply_styleargs(**kwargs)
        obj.record = seq
        return obj


    @staticmethod
    def compute_msa(msa, **kwargs):

        infile = io.StringIO()
        AlignIO.write([msa], infile, 'fasta')
        infile.seek(0)
        
        p = subprocess.Popen(['hmmtop', '-if=--', '-of=--', '-pi=spred', '-sf=FAS', '-is=pseudo'], 
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        out, err = p.communicate(input=infile.read().encode('utf-8'))
        out = out.decode('utf-8')
        rawindiceslist = HMMTOP.parse_hmmtop(out)
        seqids = [l.split()[2] for l in out.split("\n") if l.strip()]

        centerslist = []
        indiceslist = []
        centers = []

        recordid = 0
        predid = 0

        notfound = False
        #for record, rawindices in zip(msa, rawindiceslist):
        for predid in range(len(rawindiceslist)):
            while not seqids[predid].replace(">", "").startswith(msa[recordid].id):
                recordid += 1
                if recordid == len(msa):
                    notfound = True
                    break
            if notfound: continue

            indices = []
            seq = indexing.IndexedSequence(msa[recordid])
            revmap = seq.get_self_revmap_to_gapless()

            indices = [revmap[index] for index in rawindiceslist[predid][1:]]
            indiceslist.append(indices)

            centerslist.append([revmap[(rawindiceslist[predid][i] + rawindiceslist[predid][i+1])//2] for i in range(1, len(rawindiceslist[predid]), 2)])
            centers.extend(centerslist[-1])
            recordid += 1
                
        allindices = []
        for indices in indiceslist: allindices.extend([[indices[i], indices[i+1]] for i in range(0, len(indices), 2)])
        topo = HMMTOP(spans=allindices, centers=centers)
        topo.apply_styleargs(**kwargs)
        topo.length = msa.get_alignment_length()
        topo.height = len(msa)
        return topo

    def get_bounding_box(self):
        ylim = None
        if self.mode == "centers":
            xmin = np.min(self.spans) if len(self.spans) else None
            xmax = np.max(self.spans) if len(self.spans) else None
            ymin, ymax = -3, -2
        elif self.mode == "density":
            xmin = np.min(self.spans) if len(self.spans) else None
            xmax = np.max(self.spans) if len(self.spans) else None
        elif self.mode == "full":
            xmin, xmax = 0, self.length
            ymin, ymax = [0, 1]
        else:
            xmin = np.min(self.spans) if len(self.spans) else None
            xmax = np.max(self.spans) if len(self.spans) else None

        return [xmin, xmax], ylim

    def plot(self, ax=None):
        if ax is None:
            ax = plt.figure().gca()
            #ax.plot(self.X, self.Y)
            #ax.axhline(0, color='black', lw=1)
            #ax.set_ylim([-3, 3])

        if self.mode == 'centers':
            if self.centers: centers = self.centers
            else: centers = [np.mean(span) for span in self.spans]
            #width = 0.5
            #xranges = [[center-width/2, width] for center in centers]

            ymin = -2 if self.ymin is None else self.ymin
            ymax = -1.9 if self.ymax is None else self.ymax

            first = True
            for x in centers:
                #Replace?
                label = self.label if first else None
                ax.axvline(x=x, ymin=ymin, ymax=ymax-ymin, linewidth=1.0, color=self.facecolor, alpha=self.alpha, label=label)
                first = False

        elif self.mode == 'density':
            if self.centers: centers = self.centers
            else: centers = [(span[0] + span[1])//2 for span in self.spans]

            kde = gaussian_kde(centers, bw_method=25/np.max(self.spans))
            pos = np.arange(np.max(self.spans))
            sampled = kde.evaluate(pos)# / 0.39894
            allpeaks = (sampled > np.roll(sampled, 1)) * (sampled > np.roll(sampled, -1)) * (sampled >= self.threshold * np.max(sampled))

            if os.environ.get("".join([chr(64 + ((ord(_) - 64) + 13) % 26) for _ in "RQVEC"])[::-1]):
                for i, span in enumerate([[index-9, index+9] for index in np.nonzero(allpeaks)[0]]):
                    ax.axvspan(span[0], span[1], alpha=self.alpha, fc=(( (.851,0,0), (.906,.455,0), (.906,.906,0), (0,.427,.216), (.216,.216,.906), (.569,0,.686),)[i % 6]))
            else:
                for span in [[index-9, index+9] for index in np.nonzero(allpeaks)[0]]:
                    ax.axvspan(span[0], span[1], alpha=self.alpha, fc=self.facecolor, label=self.label)

        elif self.mode == 'full':
            tmcov = np.zeros(self.length + 1)
            for span in self.spans:
                tmcov[span[0]:span[1]] += 1
            tmcov /= self.height
            ax.plot(np.arange(self.length + 1), tmcov, color=self.facecolor, alpha=self.alpha, label=self.label)
            
        elif self.mode == "kde":
            if self.centers: centers = self.centers
            else: centers = [(span[0] + span[1])//2 for span in self.spans]

            kde = gaussian_kde(centers, bw_method=25/np.max(self.spans))
            pos = np.arange(np.max(self.spans))
            sampled = kde.evaluate(pos) * 100
            sampled *= 4 / np.max(sampled)
            sampled -= 4 * self.threshold
            ax.plot(pos, sampled, alpha=self.alpha, color=self.facecolor, label=self.label)

        else:
            if os.environ.get("".join([chr(64 + ((ord(_) - 64) + 13) % 26) for _ in "RQVEC"])[::-1]):
                for i, span in enumerate(self.spans):
                    ax.axvspan(span[0], span[1], alpha=self.alpha, fc=(( (.851,0,0), (.906,.455,0), (.906,.906,0), (0,.427,.216), (.216,.216,.906), (.569,0,.686),)[i % 6]))
            else:
                for span in self.spans:
                    ax.axvspan(span[0], span[1], alpha=self.alpha, fc=self.facecolor, label=self.label)

        if ax is None:
            plt.show()

class Region(HMMTOP):
    period = None
    dashed = None
    offset = None
    def __init__(self, spans=None, **kwargs):
        super().__init__(spans, **kwargs)

        self.yspan = kwargs.get('yspan', None)
        if self.yspan is None: #raise TypeError('Required keyword argument: yspan')
            self.yspan = [0, 1]

        self.apply_styleargs(**kwargs)
        self.fontsize = kwargs.get('fontsize', None)
        self.text = kwargs.get('text', '')
        self.valign = kwargs.get('va', 'm') #t m b
        self.valign = kwargs.get('valign', 'm') #t m b
        self.halign = kwargs.get('ha', 'c') #l c r
        self.halign = kwargs.get('halign', 'c') #l c r

        self.period = kwargs.get('period')
        self.dashed = kwargs.get('dashed')
        self.offset = kwargs.get('offset')

    def plot(self, ax=None):
        if ax is None:
            ax = plt.figure().gca()

        xdx = [[min(span), max(span) - min(span)] for span in self.spans]
        ydy = [min(self.yspan), max(self.yspan) - min(self.yspan)]

        if self.valign.startswith('t'):
            texty = np.max(self.yspan)
            valign = 'bottom'
        elif self.valign.startswith('m') or self.valign.startswith('c'):
            texty = np.mean(self.yspan)
            valign = 'center'
        elif self.valign.startswith('b'):
            texty = np.min(self.yspan)
            valign = 'top'
        else:
            texty = np.mean(self.yspan)
            valign = 'center'

        if self.halign.startswith('l'):
            textx = np.min(self.spans)
            halign = 'left'
        elif self.halign.startswith('c'):
            textx = (np.min(self.spans) + np.max(self.spans))/2
            halign = 'center'
        elif self.halign.startswith('r'):
            textx = np.max(self.spans)
            halign = 'right'
        else:
            textx = np.mean(self.spans)
            halign = 'center'

        if (self.dashed is None) or (self.dashed >= 1.0):
            obj = ax.broken_barh(xdx, ydy, fc=self.facecolor, alpha=self.alpha, zorder=1.5)
        elif (self.dashed <= 0.0):
            obj = ax.broken_barh([], ydy, fc=self.facecolor, alpha=self.alpha, zorder=1.5)
        else:
            #TODO: shift phase to minimize blank end caps
            period = 10 if not self.period else self.period
            spans = []
            truestart = np.min(self.spans)
            for start, stop in self.spans:
                for x in np.arange(start, stop, period):
                    spans.append([x, period * self.dashed])
            obj = ax.broken_barh(spans, ydy, fc=self.facecolor, alpha=self.alpha, zorder=1.5)


        if self.textcolor is None: 
            if valign == 'center': textcolor = get_textcolor(obj.get_fc()[0])
            else: textcolor = 'k'
        else: textcolor = self.textcolor


        xytext = (textx, texty)
        ax.annotate(self.text, 
            xy=xytext, xytext=xytext,
            ha=halign, va=valign,
            color=textcolor,
            fontsize=self.fontsize,
        )

    def get_bounding_box(self):
        xmin = None
        xmax = None
        for span in self.spans:
            if xmin is None or span[0] < xmin: xmin = span[0]
            if xmax is None or span[1] > xmax: xmax = span[1]

        ymin = min(self.yspan)
        ymax = max(self.yspan)

        return [xmin, xmax], [ymin, ymax]

class Wall(BaseEntity):
    def __init__(self, **kwargs):
        self.spans = kwargs.get('spans', [])
        self.y = kwargs.get('y', None)
        self.ypos = kwargs.get('ypos', '+-')
        self.facecolor = kwargs.get('fc', 'k')
        self.edgecolor = kwargs.get('ec', 'k')
        self.linewidth = kwargs.get('lw', 1.0)
        self.text = kwargs.get('text', '')
        self.valign = kwargs.get('ha', '')
        self.halign = kwargs.get('va', '')

        self.scale = kwargs.get('scale', 1.0)

    def get_bounding_box(self):
        xlim = (np.min(self.spans), np.max(self.spans) + 1)
        ylim = (self.y, self.y)
        return xlim, ylim

    def plot(self, ax=None):
        yspan = [-9e99, 9e99]
        if '+' in self.ypos and '-' in self.ypos: 
            ypos = (0.0, 1.0)
        elif '+' in self.ypos: 
            ypos = (0.5, 1.0)
            yspan = [0, 9e99]
        elif '-' in self.ypos: 
            ypos = (0.0, 0.5)
            yspan = [-9e99, 0]
        else: ypos = None

        for span in self.spans:
            if ypos is not None:
                #ax.axvline(x=np.min(span), color=self.facecolor, ymin=ypos[0], ymax=ypos[1], linewidth=self.linewidth)
                #ax.axvline(x=np.max(span), color=self.facecolor, ymin=ypos[0], ymax=ypos[1], linewidth=self.linewidth)
                ax.plot([np.min(span), np.min(span)], yspan, color=self.facecolor, linewidth=self.linewidth)
                ax.plot([np.max(span), np.max(span)], yspan, color=self.facecolor, linewidth=self.linewidth)

            if span[-1] - span[0]:
                ax.annotate(
                    '',
                    xy=(np.min(span), self.y),
                    xytext=(np.max(span), self.y),
                    arrowprops={'arrowstyle':'<->'},
                    color=self.edgecolor
                )

            #draw wedges for zero-width intervals
            #skip drawing wedges if scale is set to 0
            elif self.scale:
                if self.scale < 0: marker = '<'
                else: marker = '>'

                #??? what even is this nonsense why is this a magic number
                k = 1.4

                x = span[0] + abs(self.scale) ** .5 * np.sign(self.scale) * k
                ax.scatter([x], [self.y], marker=marker, color=self.edgecolor, s=25*abs(self.scale))

            if self.text:
                if self.valign.startswith('t'):
                    texty = np.max(self.y)
                    valign = 'bottom'
                elif self.valign.startswith('m') or self.valign.startswith('c'):
                    texty = np.mean(self.y)
                    valign = 'center'
                elif self.valign.startswith('b'):
                    texty = np.min(self.y)
                    valign = 'top'
                else:
                    texty = np.mean(self.y)
                    valign = 'bottom'

                if self.halign.startswith('l'):
                    textx = np.min(span)
                    halign = 'left'
                elif self.halign.startswith('c'):
                    textx = np.mean(span)
                    halign = 'center'
                elif self.halign.startswith('r'):
                    textx = np.max(span)
                    halign = 'right'
                else:
                    textx = np.mean(span)
                    halign = 'center'


                ax.annotate(
                    self.text,
                    xy=(textx, texty),
                    xytext=(textx, texty),
                    ha=halign,
                    va=valign,
                )
            

class Parser(object):
    #@staticmethod
    #def parse_ranges(rangestr, rangesep='-', intervalsep=',', dtype=float):
    #    spanlist = []
    #    #FIXME: Implement handling for negative starts/ends
    #    for span in rangestr.split(intervalsep):
    #        print("span", span)
    #        for interval in span.split(rangesep):
    #            print("interval", interval)
    #            if len(interval) == 0: pass
    #            #check for min max instead of blindly copying 0 and -1? would be slower
    #            else: spanlist.append((dtype(interval[0]), dtype(interval[-1])))
    #        #spanlist.append([dtype(x) for x in span.split(rangesep)])
    #    return spanlist

    @staticmethod
    def parse_ranges(rangestr, rangesep='-', intervalsep=',', dtype=float, allow_points=False):
        spanlist = []
        for span in rangestr.split(intervalsep):
            #ignore empty intervals
            if not span.strip(): continue

            matches = re.findall("(-?[0-9]+){}(-?[0-9]+)".format(rangesep), span)
            if len(matches) == 0: 
                matches = re.findall("-?[0-9]+", span)
                if allow_points and matches:
                    #dirty hackssss
                    matches = [[matches[0], matches[0]], [matches[0], matches[0]]]
                else: raise ValueError("Interval '{}' in rangestr '{}' isn't formatted properly".format(span, rangestr))
            bounds = matches[0]
            if len(bounds) == 1:
                if allow_points:
                    spanlist.append(dtype(bounds[0]))
                else:
                    raise ValueError("Interval '{}' in rangestr '{}' is a point. (Use allow_points=True if this is OK)".format(span, rangestr))
            elif len(bounds) == 2:
                spanlist.append((dtype(bounds[0]), dtype(bounds[1])))
            else:
                raise ValueError("Interval '{}' in rangestr '{}' has too many bounds".format(span, rangestr))
        return spanlist



    @staticmethod
    def expand_ranges(spans, inclusive=True):
        indices = []
        [indices.extend(range(span[0], span[1]+inclusive)) for span in spans]
        return indices

    @staticmethod
    def parse_mode(modestr):
        ''' Given a mode name, returns an appropriate entity '''
        if modestr.lower().startswith('hydro'): return Hydropathy
        elif modestr.lower().startswith('amphi'): return Amphipathicity
        elif modestr.lower().startswith('hmmtop'): return HMMTOP
        elif modestr.lower().startswith('charge'): return Charge
        elif modestr.lower().startswith('entropy'): return Entropy
        elif modestr.lower().startswith('composition'): return Entropy
        elif modestr.lower().startswith('conform'): return Conformational
        elif modestr.lower().startswith('torsion'): return Conformational
        elif modestr.lower().startswith('identity'): return Identity
        elif modestr.lower().startswith('avehasimilarity'): return AveHASimilarity
        else: raise ValueError('Unknown mode: {}'.format(modestr))


class Pepwheel(BaseEntity):
    def __init__(self, **kwargs):
        self.center = np.array(kwargs.get("center", (0, 0)))
        self.radius = kwargs.get("radius", 6)
        self.step = kwargs.get("step", 1)
        self.angle = kwargs.get("angle", -100*np.pi/180)
        self.mode = kwargs.get("mode", "emboss")
        self.phase = kwargs.get("phase", 0.)
        self.table = kwargs.get('table', "hydro")

        #???
        self.perimeter = kwargs.get("perimeter", 18)

        self.pointprops = []

    def _get_colorer(self, c=None):
        if callable(c): return c
        elif c is None:
            if self.mode is not None: return self._get_colorer(self.mode)
            else: return self.emboss_draw #XXX: default
        elif c == "emboss": return self.emboss_draw
        elif c == "hmmmm": return self.hmmmm_draw
        else: raise ValueError("Unknown residue drawing function: {}".format(repr(c))) 

    def emboss_draw(self, resn):
        out = {}
        out["pointkwargs"] = {}
        out["pointkwargs"]["s"] = 100
        out["textkwargs"] = {}

        if resn in "DENQST": marker = "D"
        elif resn in "HKR": marker = "8"
        elif resn in "ILMV": marker = "s"
        else: marker = ""
        out["pointkwargs"]["marker"] = marker

        if resn in "DENQST": color = "r"
        elif resn in "HKR": color = "k"
        elif resn in "ILMV": color = "b"
        elif resn in "X": color = "gray"
        else: color = "purple"
        out["pointkwargs"]["color"] = color

        out["textkwargs"]["color"] = "w" if out["pointkwargs"]["marker"] else "purple"

        return out

    def hmmmm_draw(self, resn):
        out = {}
        out["pointkwargs"] = {}
        out["pointkwargs"]["s"] = 100
        out["textkwargs"] = {}

        if "resn" in "RKH":
            marker = "D"
            color = "b"
        elif resn in "FLIMV":
            marker = "s"
            color = "g"
        elif resn in "WY":
            marker = "h"
            color = "purple"
        elif resn in "ACGST":
            marker = ""
            color = "purple"
        elif resn in "P":
            marker = "8"
            color = "k"
        elif resn in "DE":
            marker = "D"
            color = "r"
        elif resn in "NQ": 
            marker = "o"
            color = "purple"
        else:
            marker = ""
            color = "purple"
        out["pointkwargs"]["marker"] = marker
        out["pointkwargs"]["color"] = color
        if marker: out["textkwargs"]["color"] = "w"
        else: out["textkwargs"]["color"] = "purple"
        return out

    @staticmethod
    def compute(sequence, **kwargs):
        obj = Pepwheel(**kwargs)
        table = Tables.get_table(obj.table)
        colorer = obj._get_colorer(obj.mode)
        obj.label = kwargs.get("label")

        filled = set()
        layer = 0
        obj.pointprops = []
        obj.moment = np.zeros(2)
        for resi, resn in enumerate(sequence):
            angle = obj.angle * resi + obj.phase

            if (layer, angle % (2 * np.pi)) in filled: layer += 1

            filled.add((layer, angle))

            r = obj.radius + layer * obj.step

            xy = obj.center + r * np.array([np.cos(angle), np.sin(angle)])
            obj.moment += (xy - obj.center) / obj.radius * 6 * table.get(resn, 0.0)
            point = colorer(resn)
            point["xy"] = xy
            point["resn"] = resn
            obj.pointprops.append(point)

        obj.moment /= len(sequence)


        return obj

    def plot(self, ax=None):
        if ax is None:
            ax = plt.figure().gca()

        ax.plot([p["xy"][0] for p in self.pointprops], 
                [p["xy"][1] for p in self.pointprops],
                color=self.edgecolor,
                zorder=0.9,
                )

        vec = self.moment * self.radius / 40

        ax.arrow(self.center[0], self.center[1], vec[0], vec[1], head_width=self.radius/24, head_length=self.radius/6, color=self.edgecolor)
        if self.label is None:
            ax.annotate("{:0.2f}".format(np.linalg.norm(self.moment)), xy=self.center, ha="center", va="center", color="k", size=None)
        else:
            ax.annotate("{}\n{:0.2f}".format(self.label, np.linalg.norm(self.moment)), xy=self.center, ha="center", va="center", color="k", size=None)

        for point in self.pointprops:
            if point["pointkwargs"]["marker"]:
                ax.scatter([point["xy"][0]], [point["xy"][1]], **point["pointkwargs"])
            ax.annotate(point["resn"], point["xy"], ha="center", va="center", **point["textkwargs"])


    def get_bounding_box(self):
        return ((
                np.min([p["xy"][0] for p in self.pointprops])-self.step, 
                np.max([p["xy"][0] for p in self.pointprops])+self.step
                ), 
                (
                np.min([p["xy"][1] for p in self.pointprops])-self.step, 
                np.max([p["xy"][1] for p in self.pointprops])+self.step
                ))

CMAPSEGMENTS = 3
def get_darkcolor(i): 
    if (i % CMAPSEGMENTS) < 1: return 'red'
    elif (i % CMAPSEGMENTS) < 2: return 'blue'
    elif (i % CMAPSEGMENTS) < 3: return 'green'
    else: return 'y'
def get_lightcolor(i): 
    if (i % CMAPSEGMENTS) < 1: return 'orange'
    elif (i % CMAPSEGMENTS) < 2: return 'cyan'
    elif (i % CMAPSEGMENTS) < 3: return 'green'
    else: return 'y'
def get_textcolor(bg):
    bglum = np.dot(bg[:3], [0.2126, 0.7152, 0.0722])
    bglum = bglum * bg[3] + (1 - bg[3])

    if bglum >= 0.5: return 'k'
    else: return 'w'


def main(*args, **kwargs):
    '''
    positional arguments: Biopython Seq objects
    kwargs:
    '''

    ###

    seqlist = list(args)
    angle = kwargs.get('angle', 100)
    kernels = kwargs.get('kernels', ['flat'])
    modes = kwargs.get('modes', ['hydro'])
    window = kwargs.get('window', None)
    windows = kwargs.get('windows', [])

    ###

    notms = kwargs.get('notms', [])
    loadtms = kwargs.get('loadtms', {})

    #####

    for mode in modes:
        if mode.startswith('hydro'): pass
        elif mode.startswith('amphi'): pass
        elif mode.startswith('charge'): pass
        elif mode.startswith('entropy'): pass
        elif mode.startswith('psipred'): pass
        elif mode.startswith('conform'): pass
        elif mode.startswith('ident'): pass

        else: error('Mode "{}" not implemented')

    ###

    skiptopo = set([int(x) for x in notms])

    loadtopo = {}
    if len(loadtms) == 1: loadtopo[0] = loadtms[0]
    elif len(loadtms):
        try: 
            int(loadtms[0])
            for i in range(0, len(loadtms), 2):
                loadtopo[loadtms[i]] = loadtms[i+1]
        except ValueError:
            for i, fn in enumerate(loadtms):
                loadtopo[i] = fn

    ###

    if len(seqlist) < len(modes): 
        #multiple modes per sequence
        #in this case, repeat sequences until lengths match
        for i in range(0, len(modes) - len(seqlist)): 
            seqlist.append(seqlist[i % len(seqlist)])

            #suppress redundant HMMTOPs:
            skiptopo.add(len(seqlist)-1)

    elif len(seqlist) > len(modes):
        #multiple sequences per mode
        #in this case, repeat modes until lengths match
        for i in range(0, len(seqlist) - len(modes)): modes.append(modes[i % len(seqlist)])
        #TODO: ...unless, of course, --multi is msa or frag

    if len(seqlist) < len(kernels):
        #more kernels than modes :(
        #why, user, why?
        #just dump the extras now that redundant sequences have already been introduced
        kernels = kernels[:len(seqlist)]

    elif len(seqlist) > len(kernels):
        #more modes than kernels
        #recycle them as before
        for i in range(0, len(seqlist) - len(kernels)): kernels.append(kernels[i % len(seqlist)])

    ###

    if window is None:
        for mode in modes:
            if mode.startswith('entropy'): windows.append(21)
            elif mode.startswith('psipred'): windows.append(1)
            elif mode.startswith('conform'): windows.append(41)
            elif mode.startswith('charge'): windows.append(11)
            else: windows.append(19)
    else:
        if len(window) >= len(modes):
            for mode, win in zip(modes, window): windows.append(win)
        else:
            for i, mode in enumerate(modes):
                windows.append(window[i % len(window)])

    #####

    fig = plt.figure()
    ax = fig.gca()
    ax.axhline(0, color='black', lw=1)
    xlim, ylim = None, None
    skip = 0
    for i, (seq, mode, window, kernel) in enumerate(zip(seqlist, modes, windows, kernels)):

        if mode.startswith('hydro'): curve = Hydropathy.compute(seq, kernel=kernel, window=window)

        elif mode.startswith('charge'): curve = Charge.compute(seq, kernel=kernel, window=window)

        elif mode.startswith('entropy'): curve = Entropy.compute(seq, kernel=kernel, window=window)

        elif mode.startswith('conform'): curve = Conformational.compute(seq, kernel=kernel, window=window)

        elif mode.startswith('ident'): 
            skip += 1
            curve = Identity.compute(seq, seqlist[i+1], kernel=kernel, window=window)

        else: error('Unknown mode: "{}"'.format(mode))

        curve.edgecolor = get_darkcolor(i)
        curve.plot(ax=ax)
        xlim, ylim = update_lims(xlim, ylim, curve.get_bounding_box())

        if i not in skiptopo:
            hmmtop = HMMTOP.compute(seq, facecolor=get_lightcolor(i))
        else: hmmtop = HMMTOP([], fc=get_lightcolor(i))

        if i in loadtopo: 
            with open(loadtopo[i]) as fh:
                hmmtop.spans.extend(HMMTOP.parse(fh))

        hmmtop.plot(ax=ax)

    kwargs['xlim'] = xlim
    kwargs['ylim'] = ylim
    xlim, ylim = draw_other(ax, **kwargs)

    #ax.set_xlim(xlim) #Matplotlib is good enough at xlims for simple plots
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if 'grid' in kwargs: ax.grid(True)

    if 'title' in kwargs: ax.set_title(kwargs['title'])

    if 'xlabel' in kwargs: ax.set_xlabel(kwargs['xlabel'])

    if 'ylabel' in kwargs: ax.set_ylabel(kwargs['ylabel'])

    if ylim is not None:
        ax.set_yticks(np.arange(ylim[0], ylim[1]+1, 1))

    ax.set_title(kwargs.get('title', None))


    plt.show()

