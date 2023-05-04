#!/usr/bin/env python
from __future__ import print_function

import matplotlib
matplotlib.use("Agg")
import time
import os
import gzip
import sys
import argparse
import subprocess
import tempfile
import pathlib
import json
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq, Align
import Bio.Blast.NCBIXML as NCBIXML
import libquod
import fasta_grepper

try: import hvordan_gui
except ImportError: print('[WARNING]', 'Could not import hvordan_gui. Interactive interfaces will be unavailable.', file=sys.stderr)

try: from urllib.request import urlopen
except ImportError: from urllib import urlopen

def error(*things): 
    print('[ERROR]', *things, file=sys.stderr)
    exit(1)
def warn(*things): print('[WARNING]', *things, file=sys.stderr)
def info(*things): print('[INFO]', *things, file=sys.stderr)

def notfound(accfrag, acclist):
    for acc in acclist:
        if accfrag in acc: return False
        elif '-' in accfrag and '-' in acc:
            if accfrag[:accfrag.find('-')] in acc[:acc.find('-')]:
                return False
    return True

def fuzzyget(accfrag, acclist, replace=True):
    for i, acc in enumerate(acclist):
        if accfrag in acc:
            if not replace: acclist.pop(i)
            return acc
        elif '-' in accfrag and '-' in acc:
            if accfrag[:accfrag.find('-')] in acc[:acc.find('-')]:
                if not replace: acclist.pop(i)
                return acc
    return None

class IndexMapping(object):
    def __init__(self):
        self.offsets = {}

    def add(self, index, value):
        self.offsets[index] = value

    def __getitem__(self, index):
        if not self.offsets: return index
        else: return index + sum([self.offsets[_] for _ in self.offsets if _ <= index])

    def unfold(self, length):
        offset = 0
        out = {}
        for i in range(length):
            if i in self.offsets: offset += self.offsets[i]
            out[i] = i + offset
        return out

    def __call__(self, other):
        print(self.offsets, other.offsets)
        return self

    def __add__(self, other):
        out = IndexMapping()
        for index in self.offsets: out.add(index, self.offsets[index])
        for index in other.offsets:
            if index in out.offsets: out.offsets[index] += other.offsets[index]
            else: out.add(index, other.offsets[index])
        return out


class Proto1Hit(object):
    query = None
    subject = None
    bitscore = None
    expect = None
    pident = None
    qlim = None
    qlength = None
    slim = None
    slength = None

    tcblast = None
    def __init__(self, query=None, subject=None, bitscore=None, expect=None, pident=None, qlim=None, qlength=None, slim=None, slength=None):
        self.query = query
        self.subject = subject
        self.bitscore = bitscore
        self.expect = expect
        self.pident = pident
        self.qlim = qlim
        self.qlength = qlength
        self.slim = slim
        self.slength = slength

    @staticmethod
    def from_p1line(line):
        sl = line.split('\t')
        query = SeqIO.SeqRecord(Seq.Seq(sl[13].strip()), id=sl[0], name='', description='')
        subject = SeqIO.SeqRecord(Seq.Seq(sl[14].strip()), id=sl[1], name='', description='')
        bitscore = float(sl[3])
        expect = float(sl[4])
        pident = float(sl[5])
        qlim = int(sl[6]), int(sl[7])
        qlength = int(sl[8])
        slim = int(sl[9]), int(sl[10])
        slength = int(sl[11])
        return Proto1Hit(query, subject, bitscore, expect, pident, qlim, qlength, slim, slength)

    @staticmethod
    def from_p1table(fh, targets=None):
        alignments = []
        for line in fh:
            if line.startswith('#'): continue
            else:
                aln = Proto1Hit.from_p1line(line)
                if (targets is not None) and len(targets):
                    if aln.query.id in targets or aln.subject.id in targets: 
                        alignments.append(aln)
                    else: del aln
        return alignments

class Proto2Hit(object):
    subject = None
    target = None
    swzscore = None
    gsatzscore = None
    subdir = None

    ab = None
    dc = None

    blastb = None
    blastc = None

    def __init__(self, subject, target, swzscore, gsatzscore, subdir='hvordan'):
        self.subject = subject
        self.target = target
        self.swzscore = swzscore
        self.gsatzscore = gsatzscore
        self.subdir = subdir

    @staticmethod
    def from_p2line(line, subdir='hvordan'):
        sl = line.split('\t')
        subject = SeqIO.SeqRecord(sl[6], id=sl[0], name=sl[0], description=sl[0])
        target = SeqIO.SeqRecord(sl[7], id=sl[1], name=sl[1], description=sl[1])
        swzscore = float(sl[2])
        gsatzscore = float(sl[3])
        return Proto2Hit(subject, target, swzscore, gsatzscore, subdir=subdir)

    @staticmethod
    def from_p2table(fh, subdir='hvordan'):
        alignments = []
        for line in fh:
            if line.startswith('#'): continue
            elif line.startswith('Subject'): continue
            else:
                alignments.append(Proto2Hit.from_p2line(line, subdir=subdir))
        return alignments

    def write(self, fh):
        self._write_header(fh)

        self._write_plots(fh)

        self._write_footer(fh)

    def _write_header(self, fh):
        fh.write("""<!DOCTYPE html>\n<html><head><title>HVORDAN: {self.subject.id} vs {self.target.id}</title>
<link rel="stylesheet" type="text/css" href="assets/nice.css" />
<script src="assets/openclose.js"></script>
</head><body>
<h1>HVORDAN: {self.subject.id} vs {self.target.id}</h1>

<h2><button class="showhide" id="summarysh" onclick="toggle_selection('summary', 'summarysh')">-</button>
Summary</h2>
<div class="whataln" id="summary">
SS Z-score: {self.swzscore}<br/>
GSAT Z-score: {self.gsatzscore}<br/>
""".format(self=self))

    def _write_plots(self, fh):
        fh.write("""
<h2><button class="showhide" id="tcblastsh" onclick="toggle_section('tcblast', 'tcblastsh')">-</button>
TCBLAST</h2>
<div class="tcblast" id="tcblast"><a name="tcsummary"><h3>TCBLAST Summary</h3></a>

<div class="resizeable bluebars"><div class="scrollable tabular1">
<img class="bluebarplot" src="../blasts/{b}.png"/>
</div> <div class="scrollable tabular2">
<img class="bluebarplot" src="../blasts/{c}.png"/>
</div></div>

<div class="resizeable bluebars">
<div class="clear resizeable"></div><a name="pairwise"><h3>Pairwise</h3></a>
<div class="scrollable tabular1 pairwise"><hr/>
{blasttblb}
<hr/><pre>
{blastpairwiseb}
</pre>
</div>
<div class="scrollable tabular2 pairwise"><hr/>
{blasttblc}
<hr/><pre>
{blastpairwisec}
</pre>
</div>
</div>
</div>


<div class="clear"></div><a name="abcd">
<h3><button class="showhide" id="abcdsh" onclick="toggle_section('abcd', 'abcdsh')">-</button>
ABCD hydropathy plots</h3></a>
<div class="whatall" id="abcd">
    <div class="tabular1">
        A<br/><img class="bluebarplot" id="plota" src="../graphs/{aba}.png"/><br/>
        B<br/><img class="bluebarplot" id="plota" src="../graphs/{abb}.png"/><br/>
    </div>
    <div class="tabular2">
        D<br/><img class="bluebarplot" id="plota" src="../graphs/{dcd}.png"/><br/>
        C<br/><img class="bluebarplot" id="plota" src="../graphs/{dcc}.png"/><br/>
    </div>
</div>
<div class="clear"></div><br/>
<h3><button class="showhide" id="bcsh" onclick="toggle_section('bc', 'bcsh')">-</button>
<a name="bc">BC hydropathy plot</a></h3>
<div class="resizeable whataln" id="bc"><div class="scrollable"><img class="bluebarplot" id="plotbc" src="../graphs/{bc}.png"/><br/>
</div>
</div>

""".format(
    aba="{}_{}_{}".format(shorten(self.ab.query.id), shorten(self.ab.subject.id), shorten(self.ab.query.id)),
    abb="{}_{}_{}".format(shorten(self.ab.query.id), shorten(self.ab.subject.id), shorten(self.ab.subject.id)),
    dcd="{}_{}_{}".format(shorten(self.dc.query.id), shorten(self.dc.subject.id), shorten(self.dc.query.id)),
    dcc="{}_{}_{}".format(shorten(self.dc.query.id), shorten(self.dc.subject.id), shorten(self.dc.subject.id)),
    bc="{}_{}".format(shorten(self.subject.id), shorten(self.target.id)),
    b="{}".format(shorten(self.subject.id)),
    c="{}".format(shorten(self.target.id)),
    blasttblb=self.blasttbl(self.ab),
    blasttblc=self.blasttbl(self.dc),
    blastpairwiseb=self.blastpairwise(self.ab),
    blastpairwisec=self.blastpairwise(self.dc),
    ))

    def blasttbl(self, p1hit):
        out = '<table class="summtbl">\n'
        for subji, subj in enumerate(p1hit.tcblast):
            tcid = subj[:subj.find('-')]
            tccomp = subj[len(tcid)+1:]
            fam = tcid[:tcid.find('.', 2)]
            hsp = p1hit.tcblast[subj][0]

            out += '<tr class="{clsname}">'.format(clsname="oddrow" if subji % 2 else "evenrow")
            out += '<td><a href="https://www.tcdb.org/search/result.php?acc={acc}">{acc}</a></td>'.format(acc=tccomp)
            out += '<td></td>' #put TMSs here
            out += '<td><a href="https://www.tcdb.org/search/result.php?tc={fam}">{acc}</a></td>'.format(acc=tcid, fam=fam)
            out += '<td><a href="#{query}_{acc}">{hsp.bits}</a></td>'.format(query=p1hit.subject.id,
                    acc=tccomp,
                    hsp=hsp)
            out += '<td>{}</td>'.format('{:01.0e}'.format(hsp.expect)[1:])
            out +='</tr>\n'
        out += '</table>\n'
        return out

    def blastpairwise(self, p1hit):
        out =""
        for subji, subj in enumerate(p1hit.tcblast):
            tcid = subj[:subj.find('-')]
            tccomp = subj[len(tcid)+1:]
            fam = tcid[:tcid.find('.', 2)]
            hsp = p1hit.tcblast[subj][0]
            query = p1hit.subject.id

            out += """<hr/><a name="{query}_{tccomp}"></a>
>{query}
Length={qlength}

 Score = {hsp.bits:0.0f} ({hsp.score}),  Expect = {hsp.expect:01.0e}
 Identities = {hsp.identities}/{hsp.align_length} ({pident:0.0%}), Positives = {hsp.positives}/{hsp.align_length} ({psimil:0.0%}), Gaps = {hsp.gaps}/{hsp.align_length} ({pgaps:0.0%})

{hspaln}
""".format(
        query=query, tccomp=tccomp,
        qlength=len(p1hit.subject),
        hsp=hsp,
        pident=hsp.identities/hsp.align_length,
        psimil=hsp.positives/hsp.align_length,
        pgaps=hsp.gaps/hsp.align_length,
        hspaln=self.hspaln(hsp)
        )
        return out

    def hspaln(self, hsp):
        qindex = hsp.query_start
        sindex = hsp.sbjct_start
        out = ''
        for i in range(0, hsp.align_length, 60):
            qseg = hsp.query[i:i+60]
            mseg = hsp.match[i:i+60]
            sseg = hsp.sbjct[i:i+60]
            nextqindex = qindex + len(qseg) - qseg.count('-')
            nextsindex = sindex + len(sseg) - sseg.count('-')
            out += 'Query  {qindex}{qseg}  {nextqindex}\n'.format(qindex=str(qindex).ljust(5), qseg=qseg, nextqindex=nextqindex)
            out += '<span class="red">            {mseg}</span>\n'.format(mseg=mseg)
            out += 'Sbjct  {sindex}{sseg}  {nextsindex}\n'.format(sindex=str(sindex).ljust(5), sseg=sseg, nextsindex=nextsindex)
            out += '\n'
        return out

    def _write_footer(self, fh):
        fh.write("""<hr/><div class="miscinfo">Generated {timestamp}</div>
</body></html>""".format(timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ%z")))

class FullTransitivity(object):
    ab = None
    bc = None
    dc = None
    def __init__(self, ab, bc, dc):
        self.ab = ab
        self.bc = bc
        self.dc = dc

class Hvordan(object):
    p1dir = None
    p2dir = None
    outdir = None

    alignments = None

    p2alignments = None

    pfam = None

    dpi = 100

    def __init__(self, p1dir, p2dir, outdir, pfamdb):
        self.p1dir = p1dir
        self.p2dir = p2dir
        self.outdir = outdir
        self.pfamdb = pathlib.Path(pfamdb)

        self.p1alignments = []
        self.p2alignments = []
        self.full_sequences = {}

    def get_p2alignment(self, acc1, acc2=None):
        out = []
        for alignment in self.p2alignments:
            alnids = [alignment.subject.id, alignment.target.id]
            if acc2 is None:
                if not notfound(acc1, alnids): out.appent(alignment)
            else:
                if (not notfound(acc1, alnids)) and (not notfound(acc2, alnids)):
                    return alignment
        if acc2 is None: return out
        else: return None

    def get_p1alignment(self, acc1, acc2=None):
        out = []
        for alignment in self.p1alignments:
            alnids = [alignment.query.id, alignment.subject.id]
            if acc2 is None:
                if not notfound(acc1, alnids): out.append(alignment)
            else:
                if (not notfound(acc1, alnids)) and (not notfound(acc2, alnids)):
                    return alignment
        if acc2 is None: return out
        else: return None

    def collect_p2alignments(self, pairs=None, include=None, zlim=None, incfams=None, fampairs=None):
        if (pairs is not None) and (include is not None): warn('Pair list and include list are both defined! Results will be drawn from the union of these two')

        if ((pairs is not None) or (include is not None)) and (zlim is not None): warn('Accession list(s) and Z limits are both defined! Be warned that accession list(s) take precedence over Z limits!')
        alignments = []
        for path in self.p2dir:
            print(path)
            if path.name == 'report.tbl':
                with open(path) as fh: 
                    alignments.extend(Proto2Hit.from_p2table(fh))
            elif path.is_dir():
                for subpath in path.glob('*'): 
                    if incfams is not None:
                        if sum([dirfam in incfams for dirfam in subpath.name.split('_vs_')]) < 1: continue
                    if fampairs is not None:
                        if tuple(subpath.name.split('_vs_')) not in fampairs and tuple(subpath.name.split('_vs_')[::-1]) not in fampairs: continue
                    if subpath.name == 'report.tbl':
                        with open(subpath) as fh: 
                            alignments.extend(Proto2Hit.from_p2table(fh, subdir=path.name))
                    elif subpath.is_dir():
                        for subsubpath in subpath.glob('*'): 
                            if subsubpath.name == 'report.tbl':
                                with open(subsubpath) as fh: 
                                    alignments.extend(Proto2Hit.from_p2table(fh, subdir=subpath.name))
        if pairs is not None: remainingpairs = pairs[:]
        else: remainingpairs = None

        finalalignments = []
        for aln in alignments:
            found = False
            if pairs is not None:
                for i, pair in enumerate(remainingpairs):
                    if aln.subject.id in pair and aln.target.id in pair:
                        finalalignments.append(aln)
                        remainingpairs.pop(i)
                        found = True
                        break
            if found: continue

            if include is not None:
                if aln.subject.id in include or aln.target.id in include:
                    finalalignments.append(aln)
                    continue

            if zlim is not None:
                if zlim[0] is not None and aln.gsatzscore < zlim[0]: continue
                if zlim[1] is not None and aln.gsatzscore >= zlim[1]: continue

            finalalignments.append(aln)

        del alignments

        self.p2alignments = finalalignments
        return self.p2alignments

    def collect_p1alignments(self):
        targets = set()
        alignments = []
        for aln in self.p2alignments:
            targets.add(aln.subject.id)
            targets.add(aln.target.id)
        for path in self.p1dir:
            if path.name == 'psiblast.tbl':
                with open(path) as fh:
                    alignments.extend(Proto1Hit.from_p1table(fh, targets=targets))
            elif path.is_dir:
                for subpath in path.glob('*'):
                    if subpath.name == 'psiblast.tbl':
                        with open(subpath) as fh: 
                            alignments.extend(Proto1Hit.from_p1table(fh, targets=targets))
                    elif subpath.is_dir():
                        for subsubpath in subpath.glob('*'): 
                            if subsubpath.name == 'psiblast.tbl':
                                with open(subsubpath) as fh: 
                                    alignments.extend(Proto1Hit.from_p1table(fh, targets=targets))

        finalalignments = {}
        for aln in alignments:
            if aln.subject.id in finalalignments:
                if aln.expect < finalalignments[aln.subject.id].expect:
                    finalalignments[aln.subject.id] = aln
            else: finalalignments[aln.subject.id] = aln

        self.p1alignments = [finalalignments[k] for k in finalalignments]

        return self.p1alignments

    def get_full_transitivity(self):
        for alignment in self.p2alignments:
            ablist = self.get_p1alignment(alignment.subject.id)
            dclist = self.get_p1alignment(alignment.target.id)
            if len(ablist) == 1 and len(dclist) == 1:
                alignment.ab = ablist[0]
                alignment.dc = dclist[0]
            else: print("Ambiguous or missing fxp hits:", alignment.subject.id, alignment.target.id, ablist, dclist)

    def relate_p2_p1(self):
        for al2 in self.p2alignments:
            pass

    def fetch_full_sequences(self, email, force=False):
        fetchme = set()
        tcids = set()

        prefetched = set()
        (self.outdir/'_working_files').mkdir(exist_ok=True)
        if (self.outdir/'_working_files/allseqs.faa').exists():
            with open(self.outdir/'_working_files/allseqs.faa') as fh:
                prefetched = [record.id for record in SeqIO.parse(fh, 'fasta')]
        
        for aln in self.p1alignments:
            if aln.query.id.count('.') <= 1 and (force or notfound(aln.query.id, prefetched)): 
                fetchme.add(aln.query.id)
            else: tcids.add(aln.query.id)

            if force or notfound(aln.subject.id, prefetched): 
                fetchme.add(aln.subject.id)

        if force or not (self.outdir/'_working_files/tcdb.faa').exists():
            subprocess.call(['extractFamily.pl', '-i', 'all', '-o', str(self.outdir/'_working_files'), '-f', 'fasta'])

        fetchme = sorted(fetchme)

        if fetchme:
            mode = 'w' if not prefetched else 'a'
            with open(self.outdir/'_working_files/allseqs.faa', mode) as seqfh:
                for i in range(0, len(fetchme), 1000):
                    batch = fetchme[i:i+1000]
                    postdata = 'db=protein&retmode=text&rettype=fasta&email={email}&id={acclist}'.format(email=email, acclist=','.join(batch))
                    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

                    with urlopen(url, data=postdata.encode('utf-8')) as fh:
                        seqfh.write(fh.read().decode('iso-8859-1'))

                seqfh.flush()

        seqlist = list(SeqIO.parse(self.outdir/'_working_files/allseqs.faa', 'fasta'))
        acclist = [record.id for record in seqlist]
        seqdict = dict(zip(acclist, seqlist))
        tclist = list(SeqIO.parse(self.outdir/'_working_files/tcdb.faa', 'fasta'))
        tcacclist = [record.id for record in tclist]
        tcdict = dict(zip(tcacclist, tclist))
        if len(seqlist) < len(fetchme): warn('Could not find all full-length sequences! Falling back to famXpander-saved sequences')

        self.full_sequences = {}
        used_sequences = []
        for aln in self.p1alignments:
            if aln.query.id.count('.') <= 1: 
                qacc = fuzzyget(aln.query.id, acclist, replace=False)
                if qacc is not None:
                    self.full_sequences[aln.query.id] = ('entrez', qacc, seqdict[qacc])
                    used_sequences.append(seqdict[qacc])
            else: 
                qacc = fuzzyget(aln.query.id, tcacclist, replace=False)
                if qacc is not None:
                    self.full_sequences[aln.query.id] = ('tcdb', qacc, tcdict[qacc])
                    used_sequences.append(tcdict[qacc])

            sacc = fuzzyget(aln.subject.id, acclist, replace=False)
            if sacc is not None:
                self.full_sequences[aln.subject.id] = ('entrez', sacc, seqdict[sacc])
                used_sequences.append(seqdict[sacc])

        #explicitly bypasses kwarg force because it's never bad to have a fresh usedseqs.faa
        with open(self.outdir/'_working_files/usedseqs.faa', 'w') as fh:
            SeqIO.write(used_sequences, fh, 'fasta')

        return self.full_sequences

    def mkdir(self):
        self.outdir.mkdir(exist_ok=True)
        for alignment in self.p2alignments:
            (self.outdir/alignment.subdir).mkdir(exist_ok=True)
            (self.outdir/alignment.subdir/'blasts').mkdir(exist_ok=True)
            (self.outdir/alignment.subdir/'graphs').mkdir(exist_ok=True)
            (self.outdir/alignment.subdir/'html').mkdir(exist_ok=True)
            (self.outdir/alignment.subdir/'html/assets').mkdir(exist_ok=True)
            (self.outdir/alignment.subdir/'pfam').mkdir(exist_ok=True)
            (self.outdir/alignment.subdir/'sequences').mkdir(exist_ok=True)

    def run_pfam(self, force=False):
        if force or not (self.outdir/'_working_files/usedseqs.pfam').exists():
            cmd = ['hmmscan', '--cpu', '4', '--noali', '--cut_ga', '-o', '/dev/null', '--domtblout', str(self.outdir/'_working_files/usedseqs.pfam'), str(self.pfamdb), str(self.outdir/'_working_files/usedseqs.faa')]
            subprocess.call(cmd)

        #TODO: implement delta Pfams
        if force or not (self.outdir/'_working_files/pfam.json.gz').exists():
            obj = {}
            with open(self.outdir/'_working_files/usedseqs.pfam') as fh:
                for l in fh:
                    if l.startswith('#'): continue
                    elif not l.strip(): continue

                    sl = l.split()
                    hit = {
                            'pfname': sl[0],
                            'pfacc': sl[1],
                            'pflen': sl[2],
                            'qname': sl[3],
                            'qacc': sl[4],
                            'qlen': int(sl[5]),
                            'fullexpect': float(sl[6]),
                            'fullscore': float(sl[7]),
                            'fullbias': float(sl[8]),
                            'domindex': int(sl[9]),
                            'domcount': int(sl[10]),
                            'domcexpect': float(sl[11]),
                            'domiexpect': float(sl[12]),
                            'domscore': float(sl[13]),
                            'dombias': float(sl[14]),
                            'hmmcoord': [int(sl[15]), int(sl[16])],
                            'alicoord': [int(sl[17]), int(sl[18])],
                            'envcoord': [int(sl[19]), int(sl[20])],
                            'accuracy': float(sl[21]),
                            'description': ' '.join(sl[22:]).strip(),
                            }
                    if hit['qname'] in obj: obj[hit['qname']].append(hit)
                    else: obj[hit['qname']] = [hit]
            with gzip.open(self.outdir/'_working_files/pfam.json.gz', 'wt') as fh:
                self.pfam = json.dump(obj, fh, indent=1)

        with gzip.open(self.outdir/'_working_files/pfam.json.gz', 'rt') as fh:
            self.pfam = json.load(fh)

    def align_subsequences(self, force=False):
        #p1 subsequences
        pairs = {}
        for aln in self.p1alignments:
            #a-b/d-c
            pairs[aln.query.id, aln.query.id] = {
                    "sequences": (aln.query, self.full_sequences[aln.query.id][2]),
                    }
            pairs[aln.query.id, aln.subject.id] = {
                    "sequences": (aln.query, aln.subject),
                    }

        for aln in self.p2alignments:
            #b-c
            pairs[aln.subject.id, aln.target.id] = {
                    "sequences": (aln.subject, aln.target),
                    }

        for k in pairs:
            pairs[k]["mapping"] = self._get_mapping(*pairs[k]["sequences"])

    def _get_mapping(self, seq1, seq2, identical=False):
        if seq1.id == seq2.id: #please never have a non-identical sequence pair with identical accessions
            if seq1.seq == seq2.seq: return IndexMapping()
            return self._get_ident_mapping(seq1, seq2)
        else:
            aligner = Align.PairwiseAligner()
            aligner.open_gap_score = -0.5
            aligner.extend_gap_score = -0.1
            seq1gapless = Seq.Seq(str(seq1.seq).replace('-', ''))
            map1 = self._get_ident_mapping(seq1.seq, seq1gapless)
            seq2gapless = Seq.Seq(str(seq2.seq).replace('-', ''))
            map2 = self._get_ident_mapping(seq2gapless, seq2.seq)

            map12 = IndexMapping()
            for alignment in aligner.align(seq1gapless, seq2gapless):
                lastoffset = 0
                for resi1, resi2 in alignment.path:
                    if resi1 in map12.offsets: continue
                    offset = resi2 - resi1
                    if offset - lastoffset != 0: map12.add(resi1, offset - lastoffset)
                    lastoffset = offset
                break
            return map1 + map12 + map2

    def _get_ident_mapping(self, seq1, seq2):
        forward = IndexMapping()
        offset = 0
        for resi1, resn1 in enumerate(seq1):
            if resi1 + offset >= len(seq2): pass
            elif resn1 == seq2[resi1 + offset]: continue
            elif resn1 == "-": 
                forward.add(resi1, -1)
                offset -= 1
            elif seq2[resi1 + offset] == "-":
                forward.add(resi1, +1)
                offset += 1
            #else: continue (mismatches are silently accepted)
        return forward

    def run_tcblast(self, force=False):
        fxp_prots = []
        seen = set()
        for alignment in self.p1alignments:
            if alignment.subject.id in seen: continue
            fxp_prots.append(self.full_sequences[alignment.subject.id][2])
            seen.add(alignment.subject.id)
        with open(self.outdir/'_working_files/blastme.faa', 'w') as fh:
            SeqIO.write(fxp_prots, fh, 'fasta')

        if force or not (self.outdir/'_working_files/tcblast.xml').exists():
            cmd = ['blastp', '-query', str(self.outdir/'_working_files/blastme.faa'), 
                    '-subject', str(self.outdir/'_working_files/tcdb.faa'),
                    '-out', str(self.outdir/'_working_files/tcblast.xml'),
                    '-evalue', '1.0', '-gapopen', '11', '-gapextend', '1', '-matrix', 'BLOSUM62', 
                    '-comp_based_stats', 'no', '-seg', 'no', '-outfmt', '5', 
                    ]
            subprocess.call(cmd)


        p1alnlist = [aln for aln in self.p1alignments]
        self.tcblast_alignments = {}
        need_hmmtop = set()
        need_hmmtop_tcdb = set()
        with open(self.outdir/'_working_files/tcblast.xml') as fh:
            for blast in NCBIXML.parse(fh):
                self.tcblast_alignments[blast.query] = {}

                for i, aln in enumerate(p1alnlist):
                    if not notfound(aln.subject.id, [blast.query]): 
                        aln.tcblast = self.tcblast_alignments[blast.query]
                        p1alnlist.pop(i)
                        break

                need_hmmtop.add(blast.query)
                for alignment in blast.alignments:
                    if alignment.hit_id not in self.tcblast_alignments[blast.query]:
                        self.tcblast_alignments[blast.query][alignment.hit_id] = []
                    for hsp in alignment.hsps:
                        self.tcblast_alignments[blast.query][alignment.hit_id].append(hsp)
                    need_hmmtop_tcdb.add(alignment.hit_id)

            
        for tcid in need_hmmtop_tcdb:
            self.full_sequences[tcid] =  ("tcdb", tcid, fasta_grepper.main(tcid, open(self.outdir/'_working_files/tcdb.faa'), single=True)[0])

        hmmtop_tcdb = {}
        if force or not (self.outdir/'_working_files/fullseq_hmmtop.json.gz').exists():
            #need_hmmtop_sequences = [fasta_grepper.main(tcid, open(self.outdir/'_working_files/tcdb.faa'), single=True)[0] for tcid in need_hmmtop_tcdb]
            for tcid in need_hmmtop_tcdb:
                hmmtop_tcdb[tcid] = libquod.entities.HMMTOP.compute(self.full_sequences[tcid][2])

            hmmtop_all = {}
            for acc in self.full_sequences:
                hmmtop_all[acc] = libquod.entities.HMMTOP.compute(self.full_sequences[acc][2])
            hmmtop_all.update(hmmtop_tcdb)

            out = {}
            for acc in hmmtop_all:
                out[acc] = hmmtop_all[acc].spans
            with gzip.open(self.outdir/'_working_files/fullseq_hmmtop.json.gz', 'wt') as fh:
                json.dump(out, fh)

        with gzip.open(self.outdir/'_working_files/fullseq_hmmtop.json.gz', 'rt') as fh:
            self.fullseq_hmmtop = json.load(fh)


    def run_quod(self, force=False):
        full_sequence_plots = {}
        width = 30/2
        height = 5.5/2
        for alignment in self.p2alignments:
            if alignment.ab is None or alignment.dc is None:
                print("Missing AB/DC information for", alignment.subject.id, alignment.target.id)
            else:
                outdir = self.outdir/alignment.subdir
                #A: full sequence for A
                fig, ax = plt.subplots()
                fig.set_tight_layout(True)
                fullseq = self.full_sequences[alignment.ab.query.id]
                libquod.entities.Hydropathy.compute(fullseq[2], edgecolor="brown").plot(ax=ax)
                libquod.entities.HMMTOP.compute(fullseq[2], facecolor="darkorange").plot(ax=ax)
                libquod.entities.Wall(spans=[alignment.ab.qlim], y=-2, ypos='+-').plot(ax=ax)
                pfacc = fuzzyget(fullseq[2].id, self.pfam)
                if pfacc is not None:
                    for domain in self.pfam[pfacc]:
                        libquod.entities.Region(spans=[domain["envcoord"]], 
                                text="{pfname} ({pfacc})".format(**domain),
                                valign="t",
                                halign="c",
                                fc="red",
                                alpha=1.0,
                                yspan=[-2.8,-2.6],
                                ).plot(ax=ax)

                ax.axhline(0, lw=0.5, color='k')
                ax.set_ylim([-3, 3])
                ax.set_xlim([0, len(fullseq[2])])
                ax.set_title(alignment.ab.query.id)
                fig.set_figwidth(width)
                fig.set_figheight(height)
                fig.savefig(outdir/"graphs/{}_{}_{}.png".format(shorten(alignment.ab.query.id), shorten(alignment.ab.subject.id), shorten(alignment.ab.query.id)), dpi=self.dpi)
                plt.close()

                #B: full sequence for B
                fig, ax = plt.subplots()
                fig.set_tight_layout(True)
                fullseq = self.full_sequences[alignment.subject.id]
                libquod.entities.Hydropathy.compute(fullseq[2], edgecolor="r").plot(ax=ax)
                libquod.entities.HMMTOP.compute(fullseq[2], facecolor="orange").plot(ax=ax)
                libquod.entities.Wall(spans=[alignment.ab.slim], y=+2, ypos='+').plot(ax=ax)
                pfacc = fuzzyget(fullseq[2].id, self.pfam)
                if pfacc is not None:
                    for domain in self.pfam[pfacc]:
                        libquod.entities.Region(spans=[domain["envcoord"]], 
                                text="{pfname} ({pfacc})".format(**domain),
                                valign="t",
                                halign="c",
                                fc="blue",
                                alpha=1.0,
                                yspan=[-2.8,-2.6],
                                ).plot(ax=ax)
                ax.axhline(0, lw=0.5, color='k')
                ax.set_ylim([-3, 3])
                ax.set_xlim([0, len(fullseq[2])])
                ax.set_title(alignment.subject.id)
                fig.set_figwidth(width)
                fig.set_figheight(height)
                fig.savefig(outdir/"graphs/{}_{}_{}.png".format(shorten(alignment.ab.query.id), shorten(alignment.ab.subject.id), shorten(alignment.subject.id)), dpi=self.dpi)
                plt.close()

                #D: full sequence for D
                fig, ax = plt.subplots()
                fig.set_tight_layout(True)
                fullseq = self.full_sequences[alignment.dc.query.id]
                libquod.entities.Hydropathy.compute(fullseq[2], edgecolor="darkblue").plot(ax=ax)
                libquod.entities.HMMTOP.compute(fullseq[2], facecolor="darkcyan").plot(ax=ax)
                libquod.entities.Wall(spans=[alignment.dc.qlim], y=-2, ypos='+-').plot(ax=ax)
                pfacc = fuzzyget(fullseq[2].id, self.pfam)
                if pfacc is not None:
                    for domain in self.pfam[pfacc]:
                        libquod.entities.Region(spans=[domain["envcoord"]], 
                                text="{pfname} ({pfacc})".format(**domain),
                                valign="t",
                                halign="c",
                                fc="blue",
                                alpha=1.0,
                                yspan=[-2.8,-2.6],
                                ).plot(ax=ax)
                ax.axhline(0, lw=0.5, color='k')
                ax.set_ylim([-3, 3])
                ax.set_xlim([0, len(fullseq[2])])
                ax.set_title(alignment.dc.query.id)
                fig.set_figwidth(width)
                fig.set_figheight(height)
                fig.savefig(outdir/"graphs/{}_{}_{}.png".format(shorten(alignment.dc.query.id), shorten(alignment.dc.subject.id), shorten(alignment.dc.query.id)), dpi=self.dpi)
                plt.close()

                #C: full sequence for C
                fig, ax = plt.subplots()
                fig.set_tight_layout(True)
                fullseq = self.full_sequences[alignment.target.id]
                libquod.entities.Hydropathy.compute(fullseq[2], edgecolor="b").plot(ax=ax)
                libquod.entities.HMMTOP.compute(fullseq[2], facecolor="cyan").plot(ax=ax)
                libquod.entities.Wall(spans=[alignment.dc.slim], y=+2, ypos='+').plot(ax=ax)
                pfacc = fuzzyget(fullseq[2].id, self.pfam)
                if pfacc is not None:
                    for domain in self.pfam[pfacc]:
                        libquod.entities.Region(spans=[domain["envcoord"]], 
                                text="{pfname} ({pfacc})".format(**domain),
                                valign="t",
                                halign="c",
                                fc="blue",
                                alpha=1.0,
                                yspan=[-2.8,-2.6],
                                ).plot(ax=ax)
                ax.axhline(0, lw=0.5, color='k')
                ax.set_ylim([-3, 3])
                ax.set_xlim([0, len(fullseq[2])])
                ax.set_title(alignment.target.id)
                fig.set_figwidth(width)
                fig.set_figheight(height)
                fig.savefig(outdir/"graphs/{}_{}_{}.png".format(shorten(alignment.dc.query.id), shorten(alignment.dc.subject.id), shorten(alignment.target.id)), dpi=self.dpi)
                plt.close()

                #BC:
                fig, ax = plt.subplots()
                fig.set_tight_layout(True)
                fullseqb = self.full_sequences[alignment.subject.id]
                fullseqc = self.full_sequences[alignment.target.id]
                libquod.entities.Hydropathy.compute(seq=fullseqb[2], fragment=alignment.subject, edgecolor="r").plot(ax=ax)
                libquod.entities.HMMTOP.compute(seq=fullseqb[2], fragment=alignment.subject, facecolor="orange").plot(ax=ax)
                libquod.entities.Hydropathy.compute(seq=fullseqc[2], fragment=alignment.target, edgecolor="b").plot(ax=ax)
                libquod.entities.HMMTOP.compute(seq=fullseqc[2], fragment=alignment.target, facecolor="c").plot(ax=ax)
                ax.axhline(0, lw=0.5, color='k')
                ax.set_ylim([-3, 3])
                ax.set_xlim([0, max(len(alignment.subject), len(alignment.target))])
                ax.set_title("{} vs. {}".format(alignment.subject.id, alignment.target.id))
                fig.set_figwidth(width*2)
                fig.set_figheight(height)
                fig.savefig(outdir/"graphs/{}_{}.png".format(shorten(alignment.subject.id), shorten(alignment.target.id), dpi=self.dpi))
                plt.close()

    def plot_tcblast(self):

        def expect2alpha(expect):
            #>=1e+1: 0%
            #1e-1: 20%
            #1e-3: 40%
            #1e-5: 60%
            #1e-7: 80%
            #<=1e-8: 90%
            return np.clip(-np.log10(expect)/10 + .1, 0.01, 1.0) if expect else 1.0

        done = set()
        tcblasts = {}
        for alignment in self.p2alignments:
            for prot in alignment.subject, alignment.target:
                if prot.id in self.tcblast_alignments: 
                    k = prot.id
                else: 
                    k = fuzzyget(prot.id, self.tcblast_alignments)
                    if k is None: warn('Could not find', prot.id, 'in TCBLAST results')
                if prot.id in self.full_sequences:
                    full = self.full_sequences[prot.id][2]
                else:
                    full = fuzzyget(prot.id, self.full_sequences)[2]
                    if full is None: warn('Could not find', prot.id, 'in full sequences')

                fig, ax = plt.subplots()
                fig.set_tight_layout(True)
                fig.set_figwidth(6)
                covheight = len(self.tcblast_alignments[k])/12.0 * 3
                topheight = 0.25
                bottomheight = 1.0
                fig.set_figheight(covheight + topheight + bottomheight)
                ax.get_yaxis().set_visible(0)
                ax.xaxis.tick_top()
                ax.broken_barh([[start, end - start] for start, end in self.fullseq_hmmtop[prot.id]],
                        [0, 1],
                        fc='k')
                y = 0
                for subject in self.tcblast_alignments[k]:
                    spans = []
                    for hsp in self.tcblast_alignments[k][subject]:
                        spans.extend(self._transform_tmss(
                                self.full_sequences[subject][2],
                                prot,
                                self.fullseq_hmmtop[subject], 
                                hsp,
                        ))
                    y -= 2
                    alpha = expect2alpha(self.tcblast_alignments[k][subject][0].expect)
                    ax.broken_barh([[hsp.query_start, hsp.query_end - hsp.query_start] for hsp in self.tcblast_alignments[k][subject]], [y, 1], fc='#0099ff', alpha=alpha)
                    ax.broken_barh([[start, end - start] for start, end in spans], [y, 1], fc='#333333', alpha=alpha)
                    ax.annotate(shorten(subject), (self.tcblast_alignments[k][subject][0].query_start, y), color='w', size=8)

                ax.set_xlim([0, len(full)])
                fig.savefig(self.outdir/alignment.subdir/'blasts/{}.png'.format(shorten(prot.id)), dpi=self.dpi)
                plt.close()

    def _transform_tmss(self, full1, full2, spans, hsp):
        resi = hsp.sbjct_start

        frag1 = SeqIO.SeqRecord(seq=hsp.sbjct, id='frag1', name='', description='')
        frag2 = SeqIO.SeqRecord(seq=hsp.query, id='frag2', name='', description='')

        #get full sbjct -> frag sbjct mapping
        m1 = self._get_ident_mapping(full1, frag1)
        #get frag sbjct -> frag query mapping
        m2 = self._get_ident_mapping(frag1, frag2)
        #get frag query -> full query mapping
        m3 = self._get_ident_mapping(frag2, full2)

        mnet = m1 + m2 + m3
        out = []
        for start, end in spans:
            tmstart = np.clip(mnet[start], hsp.query_start, hsp.query_end)
            tmend = np.clip(mnet[end], hsp.query_start, hsp.query_end)
            if tmend - tmstart > 0:
                out.append([tmstart, tmend])

        return out


    def write_html(self):
        subdirs = set()
        for alignment in self.p2alignments:
            subdirs.add(alignment.subdir)
            fn = self.outdir/alignment.subdir/"html/{}_vs_{}.html".format(
                    shorten(alignment.subject.id), 
                    shorten(alignment.target.id),
                    )

            with open(fn, "w") as fh:
                alignment.write(fh)

        for subdir in subdirs:
            with open(self.outdir/subdir/"html/assets/openclose.js", "w") as fh:
                fh.write("""
function toggle_section(sectionid, selfid) {
	var section = document.getElementById(sectionid);
	var me = document.getElementById(selfid);
	//console.log([section, section.style.display]);
	if (section.style.display == 'none') {
		section.style.display = 'block';
		me.innerHTML = '-';
	} else { 
		section.style.display = 'none'; 
		me.innerHTML = '+';
	}
}""")
            with open(self.outdir/subdir/"html/assets/nice.css", "w") as fh:
                fh.write("""body {

	font-family: sans-serif;
	height: 100%;
}
div {
	display: block;
}
div.tcblast {
	max-width: 1500px;
}
div.fullblast {
	width: 50%;
	float: left;
}
div.tabular1 {
	width: 49%;
	float: left;
	height: auto;
}
div.tabular2 {
	width: 49%;
	float: right;
	height: auto; 
}
img.bluebarplot {
	max-width: 99%;
	height: auto;
}
.clear { clear: both; }
.scrollable {
	overflow-y: scroll;
}
.resizeable {
	/* resize: vertical; */
	overflow: auto;
	border: 1px solid gray;
	display: block;
	padding-bottom: 1ex;
}
.bluebars {
	height: 25vh;
}
.pairwise {
	height: 50vh;
}
.whatall {
	margin: 8pt;
	/*height: 50vh;*/
}
.whataln {
	width: 100%;
}
#seqs {
	display: none;
}
.summtbl {
	font-family: monospace, courier;
	font-size: 75%;
}
.oddrow {
	background-color: #d8d8d8;
}
td {
	padding-right: 1em;
}
.red {
	color: red;
}
img {
	border: 1pt solid black;
}
.monospace {
	font-family: monospace;
}
.miscinfo {
	font-size: 4pt;
}
""")
        


def shorten(accession):
    out = accession
    if ' ' in out: out = out[:out.find(' ')]
    if out.count('.') == 1: out = out[:out.find('.')]
    return out



def run_hvordan(args):
    if args.p:
        if len(args.p) % 2: error('Number of accessions given to -p must be even')
    if args.interactive:
        zselector = hvordan_gui.HvordanZSelector(args.p2d)
        zselector.run()
        zlim = zselector.get_zlim()
    else: zlim = (args.z_min, args.z_max)

    hvor = Hvordan(args.p1d, args.p2d, args.o, pfamdb=args.pfamdb)

    if args.v >= 1: info("Collecting protocol2 alignments")
    hvor.collect_p2alignments(pairs=args.p, include=args.i, zlim=zlim, incfams=args.include_fams, fampairs=[tuple(args.include_fam_pairs[i:i+2]) for i in range(0, len(args.include_fam_pairs), 2)] if args.include_fam_pairs is not None else None)
    if args.v >= 1: info("Collecting famXpander alignments")
    hvor.collect_p1alignments()
    hvor.relate_p2_p1()

    hvor.mkdir()
    if args.v >= 1: info("Fetching full sequences")
    hvor.fetch_full_sequences(args.email, force=args.force)

    if args.v >= 1: info("Running hmmscan")
    hvor.run_pfam()

    if args.v >= 1: info("Aligning subsequences to sequences")
    hvor.align_subsequences()

    if args.v >= 1: info("Getting full transitivity")
    hvor.get_full_transitivity()

    if args.v >= 1: info("Running TC-BLASTs")
    hvor.run_tcblast()

    if args.v >= 1: info("Generating plots")
    hvor.run_quod()
    hvor.plot_tcblast()

    if args.v >= 1: info("Writing HTML")
    hvor.write_html()
    
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--p1d', nargs='+', type=pathlib.Path, required=True, help='FamXpander tables/directory/ies')

    parser.add_argument('--p2d', nargs='+', type=pathlib.Path, required=True, help='Protocol2 tables/directory/ies')
    parser.add_argument('-o', required=True, type=pathlib.Path, help='Where to write HVORDAN summaries')

    parser.add_argument('--include-fams', nargs='+', help='Which families to require. Defaults to all')
    parser.add_argument('--include-fam-pairs', nargs='+', help='Which family pairs to require. Defaults to all')

    parser.add_argument('--z-min', type=float, default=14., help='Minimum Z-score (default: 14)')
    parser.add_argument('--z-max', type=float, default=None, help='Maximum Z-score (default: none)')

    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing output if present')

    parser.add_argument('-r', '--dpi', type=int, default=100, help='Resolution of plots (default: 100)')

    parser.add_argument('-m', '--max-hits', type=int, default=10, help='Maximum target sequences for TCBLAST')

    if 'ENTREZ_EMAIL' in os.environ:
        parser.add_argument('-e', '--email', default=os.environ['ENTREZ_EMAIL'], help='Email to provide to NCBI for Entrez requests. Defaults to $ENTREZ_EMAIL if set. (default: $ENTREZ_EMAIL == {}'.format(os.environ['ENTREZ_EMAIL']))
    else:
        parser.add_argument('-e', '--email', required=True, help='Email to provide to NCBI for Entrez requests. Defaults to $ENTREZ_EMAIL if set. (required)')

    if 'PFAMDB' in os.environ:
        parser.add_argument('-d', '--pfamdb', default=os.environ['PFAMDB'], type=pathlib.Path, help='Which profile database to use. Defaults to $PFAMDB if set. (default: $PFAMDB == {})'.format(os.environ['PFAMDB']))
    else:
        parser.add_argument('-d', '--pfamdb', required=True, type=pathlib.Path, help='Which profile database to use. Defaults to $PFAMDB if set. (required)')

    parser.add_argument('--kernel', default='flat', help='Use a kernel other than the default moving average (e.g. hann, triang, gauss)')

    parser.add_argument('-i', nargs='+', help='Operate only on pairs containing these accessions')

    parser.add_argument('-p', nargs='+', help='Operate only on these specific pairs')

    parser.add_argument('--interactive', action='store_true', help='Use interactive selector')

    parser.add_argument('-v', action='store_true', help='Verbose output')

    args = parser.parse_args()

    run_hvordan(args)

if __name__ == '__main__':
    main()
