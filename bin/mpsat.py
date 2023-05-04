#!/usr/bin/env python

from __future__ import print_function, division, unicode_literals

import os
import argparse
import socket
import time
import subprocess
import shlex
import io
import sys
import tempfile
import re
import numpy as np

import scipy.stats

from Bio import SeqIO, Seq
import shuffle_sequences_within_class as shuffleseq

def error(*things):
    print('[ERROR]', *things, file=sys.stderr)
    exit(1)

NUMBINS = 20

def get_sacc(s):
    for l in s.split('\n'):
        if l.startswith('>>>'): pass
        elif l.startswith('>>') and not l.startswith('>>shuf_'):
            return l[2:]
    raise IndexError

def parse_toplike(s):
    indices = []
    topo = 1
    for rangestr in s.split(','):
        if '-' in rangestr:
            points = [int(x) for x in rangestr.split('-')]
            indices.append([points[0], points[-1]])
        elif rangestr.lower().startswith('i'): topo = 1
        elif rangestr.lower().startswith('o'): topo = 0
    return topo, indices


def main(queryfn, subjectfn, count=1000, outdir=None, flags=None, split_inout=False, split_looptail=False, shufbuffer=1000, interactive=False, expectsmooth=5, quiet=False, binsize=None, bandwidth=None, preview=False, resume=False, compact=False, stms=None, gsat=False, pad_unalignable=False, save=-1, dpi=None, width=None, height=None):
    mode = 'swzscore'
    flags = '' if flags is None else flags

    if resume is not False:
        path = outdir if resume is None else resume

        if os.path.isdir(path): 
            outfile = None
            outdir = path
        elif os.path.isfile(path): 
            outfile = path
            outdir = None
        resume = True

    topo = None
    indices = None
    if gsat:
        topo = 1
        indices = []
    elif stms:
        topo, indices = parse_toplike(','.join(stms))

    args = []
    ignore = True
    for flag in shlex.split(flags):
        if flag.startswith('-'):
            ignore = False
            if flag.startswith('-E'): ignore = True
            elif flag.startswith('-m'): ignore = True
            elif flag.startswith('-k'): ignore = True
            if not ignore: args.append(flag)
        elif ignore: continue
        else:
            args.append(flag)
    args.extend(['-E', str(count), '-k', '1000', '-m', '10'])

    if resume:
        if outdir:
            with open('{}/results.ssearch'.format(outdir)) as f: out = f.read()
            seqlist = list(SeqIO.parse('{}/subject.faa'.format(outdir), 'fasta'))
            sacc = get_sacc(out) if not seqlist else seqlist[0].name

        elif outfile:
            with open('{}'.format(outfile)) as f: out = f.read()
            sacc = get_sacc(out)

    else:
        if outdir: 
            if not os.path.isdir(outdir): os.mkdir(outdir)
            qfh = open('{}/query.faa'.format(outdir), 'w')
            rsfh = open('{}/subject.faa'.format(outdir), 'w+')
            ssfh = open('{}/subject_shuffled.faa'.format(outdir), 'w')
        else:
            qfh = tempfile.NamedTemporaryFile('w')
            rsfh = tempfile.NamedTemporaryFile('w+')
            ssfh = tempfile.NamedTemporaryFile('w')

        if queryfn.startswith('asis:'): qfh.write(format(SeqIO.SeqRecord(seq=Seq.Seq(queryfn[5:]), id='query', name='', description=''), 'fasta'))
        else: 
            with open(queryfn) as f: qfh.write(f.read())
        if subjectfn.startswith('asis:'): qfh.write(format(SeqIO.SeqRecord(seq=Seq.Seq(subjectfn[5:]), id='subject', name='', description=''), 'fasta'))
        else:
            with open(subjectfn) as f: rsfh.write(f.read())

        qfh.flush()
        rsfh.flush()
        rsfh.seek(0)



        #not sure if this should be scaled up
        sseq = SeqIO.read(rsfh, 'fasta')
        sacc = sseq.name
        rsfh.seek(0)

        SeqIO.write([sseq], ssfh, 'fasta')
        mergedict = shuffleseq.get_mergedict(split_inout, split_looptail)

        sequences = shuffleseq.shuffle_seq(sseq, count=count, prefix='shuf', mergedict=mergedict, topo=topo, indices=indices)
        SeqIO.write(sequences, ssfh, 'fasta')
        ssfh.flush()

        out = subprocess.check_output(['ssearch36'] + args + [qfh.name, ssfh.name])
        out = out.decode('utf-8')

        if outdir: 
            with open('{}/results.ssearch'.format(outdir), 'w') as f: f.write(out)

            if indices is not None:

                with open('{}/subject.top'.format(outdir), 'w') as f:
                    f.write('>HA: {slen} {sacc}   {topo}   {indices} \n'.format(
                        slen=len(sseq),
                        sacc=sacc,
                        topo='IN' if topo == 1 else 'OUT',
                        indices='  '.join([' '.join([str(pt) for pt in span]) for span in indices])
                    ))

    stats, plotdata = get_statistic(out, accession=sacc, mode=mode, expectsmooth=expectsmooth, binsize=binsize, bandwidth=bandwidth, preview=preview or save, pad_min=pad_unalignable, shuffles=count)
    stats['mpsat'] = {
            'query': queryfn,
            'subject': subjectfn,
    }

    if not compact: fmtout = get_formatted(stats)
    else: fmtout = get_compact_formatted(stats)


    if outdir: 
        if resume and os.path.exists('{}/out.mpsat'.format(outdir)): pass
        else:
            with open('{}/out.mpsat'.format(outdir), 'w') as f: f.write(fmtout)

    if interactive: raw_input(qfh.name + ' ' + ssfh.name)
    if not quiet: print(fmtout)
    if save == -1: save = False
    elif save is None: save = '{}/plot.png'.format(outdir)

    if preview or save:
        import matplotlib
        import matplotlib.pyplot as plt
        fig = plt.figure()
        fig.set_tight_layout(True)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        for k in 'KDE', 'GEV s unshuf':
            ax1.plot(*plotdata[k], label=k.replace(' s unshuf', ''))
        ax1.bar(*plotdata['Histogram centers'], color=matplotlib.cm.tab20(0.051), label='Histogram')
        ax1.legend()
        ax1.scatter(*plotdata['Unshuffled'], zorder=2.1)

        for k in 'cov KDE', 'cov GEV s unshuf':
            ax2.plot(*plotdata[k], label=k.replace('cov ', '').replace(' s unshuf', ''))
        ax2.bar(*plotdata['cov Histogram centers'], color=matplotlib.cm.tab20(0.051), label='Histogram')
        ax2.legend()
        ax2.scatter(*plotdata['cov Unshuffled'], zorder=2.1)

        xlim = [min(ax1.get_xlim()[0], ax2.get_xlim()[0]), max(ax1.get_xlim()[1], ax2.get_xlim()[1])]
        ylim = [0, max(ax1.get_ylim()[1], ax2.get_ylim()[1])]
        ax1.set_xlim(xlim)
        ax2.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax2.set_ylim(ylim)

        if width is not None: fig.set_figwidth(width)
        if height is not None: fig.set_figheight(height)
        if save: fig.savefig(save, dpi=dpi)
        if preview: plt.show()

    return fmtout

def get_kv(line):
    k = line[2:line.find(':')]
    v = line[line.find(':')+1:].strip()

    try: v = int(v)
    except ValueError:
        try: v = float(v)
        except ValueError: pass
    return k, v


def get_statistic(text, accession, mode='swzscore', expectsmooth=5, binsize=None, bandwidth=None, preview=False, pad_min=False, shuffles=0):
    section = ''
    subsection = ''

    rank = 1
    hits = 0

    currentacc = ''
    runtime = {}
    scores = []
    bitscores = []
    special = {}
    statistics = {}
    sequence = {}
    sequence['starts'] = []
    sequence['stops'] = []
    sequence['lengths'] = []

    bestexpect = []
    allexpect = []

    plotdata = {}
                
    multiline = False
    mlkey = ''

    tmpstart = 0
    alnlengths = []

    for l in text.split('\n'):
        if l.startswith('>>>'): section = 'runtime'
        elif not section: continue
        elif l.startswith('>>'): 
            section = 'alignment'
            currentacc = l[2:].strip()
            special[currentacc] = {}
            if currentacc == accession or (not currentacc.startswith('shuf_')): statistics['rank'] = rank
            rank += 1
            multiline = False
        elif l.startswith('>'): 
            section = 'sequence'
            if currentacc.startswith(accession) or accession.startswith(currentacc) or (not currentacc.startswith('shuf_')): 
                if '_order' not in sequence: sequence['_order'] = []
                mlkey = l[1:].strip()
                sequence['_order'].append(mlkey)
                multiline = True

        #elif multiline and not l.startswith('; ') and (currentacc == accession):
        elif multiline and not l.startswith('; '):
            if mlkey in sequence: sequence[mlkey] += l.replace('\n', '')
            else: sequence[mlkey] = l.replace('\n', '')

        elif subsection: pass

        elif l.startswith(';'): 
            if section == 'runtime':
                k, v = get_kv(l)
                runtime[k] = v
            elif section == 'alignment':
                k, v = get_kv(l)
                if currentacc == accession or not currentacc.startswith('shuf_'): 
                    special[currentacc][k] = v
                elif ' opt' in k: scores.append(v)
                elif 'bits' in k: bitscores.append(v)

                if currentacc.startswith('shuf_') and 'expect' in k:
                    if (len(bestexpect) < expectsmooth) and v > 0: bestexpect.append(v)
                    allexpect.append(v)
            elif section == 'sequence':
                k, v = get_kv(l)
                if (currentacc == accession) or (not currentacc.startswith('shuf_')): 
                    if 'al_display_start' in k: sequence[mlkey] = ''
                    elif 'al_cons' in k: 
                        mlkey = 'cons'
                        sequence['_order'].append(mlkey)
                        sequence[mlkey] = ''
                    elif 'al_start' in k: sequence['starts'].append(v)
                    elif 'al_stop' in k: sequence['stops'].append(v)
                    elif 'sq_len' in k: sequence['lengths'].append(v)
                    else: special[currentacc][k] = v
                else:
                    if 'al_start' in k: tmpstart = v
                    elif 'al_stop' in k: 
                        length = v - tmpstart + 1
                        if alnlengths:
                            if alnlengths[-1][1] is None: alnlengths[-1][1] = length
                            else: alnlengths.append([length, None])
                        else:
                            alnlengths.append([length, None])

    coverages = [max([(qaligned+1)/sequence['lengths'][0], (saligned+1)/sequence['lengths'][1]]) for qaligned, saligned in alnlengths]
    covbitscores = [bits*cov for (bits, cov) in zip(bitscores, coverages)]
    covscores = [score*cov for (score, cov) in zip(scores, coverages)]

    statistics['npad'] = 0
    if pad_min and shuffles > 0:
        npad = shuffles - len(scores)
        statistics['npad'] = npad
        scores.extend([scores[-1]] * npad)
        bitscores.extend([bitscores[-1]] * npad)
        allexpect.extend([allexpect[-1]] * npad)


    if accession not in special: 
        if special:
            #FIXME: Better handling for complex FASTA headers
            k = sorted(special)[0]
            special[accession] = special[k]
        else: raise Exception('Could not find subject sequence in {} alignments'.format(len(bitscores)))


    for k in special[accession]:
        if 'opt' in k: 
            statistics['swmean'] = np.mean(scores)
            statistics['swstd'] = np.std(scores)
            statistics['swzscore'] = (special[accession][k] - statistics['swmean']) / statistics['swstd']
            statistics['count'] = len(scores)
        elif 'bits' in k: pass
        if 'opt' in k:
            #statistics['bitsmean'] = np.mean(bitscores)
            #statistics['bitsstd'] = np.std(bitscores)
            statistics['bitsmean'] = np.mean(scores)
            statistics['bitsstd'] = np.std(scores)
            statistics['bitszscore'] = (special[accession][k] - statistics['bitsmean']) / statistics['bitsstd']
            statistics['count'] = len(scores)
            statistics['coverage'] = max([(sequence['stops'][_] - sequence['starts'][_]) / sequence['lengths'][_] for _ in range(2)])

            #allbitscores = np.concatenate([bitscores, [special[accession][k]]])
            allbitscores = np.concatenate([scores, [special[accession][k]]])
            allscores = np.concatenate([scores, [special[accession][k]]])
            kde = scipy.stats.gaussian_kde(allbitscores, bw_method=bandwidth)
            #covkde = scipy.stats.gaussian_kde(covbitscores, bw_method=bandwidth)
            covkde = scipy.stats.gaussian_kde(covscores, bw_method=bandwidth)
            statistics['p-value'] = kde.integrate_box_1d(special[accession][k], 9e99)
            statistics['shortcov-p-value'] = covkde.integrate_box_1d(special[accession][k] * statistics['coverage'], 9e99)
            #REMOVED FOR OPTIMIZATION
            #statistics['p-value-2x'] = kde.integrate_box_1d(special[accession][k], 2*np.max(allbitscores))

            #histcounts, histbins = np.histogram(allbitscores, bins=NUMBINS) 
            #curarea = np.sum(histcounts * (histbins[1] - histbins[0]))
            #normcounts = histcounts / curarea

            #histX = (histbins[:-1] + histbins[1:]) / 2
            #histY = normcounts
            #histYkde = kde.evaluate(histX)

            #histkde_mae = np.mean(np.abs(histYkde - histY))
            #histkde_mse = np.mean((histYkde - histY)**2)
            #statistics['histkde_mae'] = histkde_mae
            #statistics['histkde_mse'] = histkde_mse

            rawgevparams = scipy.stats.genextreme.fit(allbitscores)
            #rawshortgevparams = scipy.stats.genextreme.fit(bitscores)
            rawshortgevparams = scipy.stats.genextreme.fit(scores)
            #rawshortcovgevparams = scipy.stats.genextreme.fit(covbitscores)
            rawshortcovgevparams = scipy.stats.genextreme.fit(covscores)
            statistics['gev-p-value'] = 1 - scipy.stats.genextreme.cdf(special[accession][k], *rawgevparams) 
            statistics['shortgev-p-value'] = 1 - scipy.stats.genextreme.cdf(special[accession][k], *rawshortgevparams) 

            #rawshortcovgevparams = scipy.stats.genextreme.fit(covbitscores)
            rawshortcovgevparams = scipy.stats.genextreme.fit(covscores)
            statistics['shortcovgev-p-value'] = 1 - scipy.stats.genextreme.cdf(special[accession][k] * statistics['coverage'], *rawshortcovgevparams) 

            a = b = 1
            #n = len(bitscores)
            n = len(scores)
            r = np.clip(statistics['rank'] - 1, 0, n)
            #rcov = sum(np.array(covbitscores) > (special[accession][k] * statistics['coverage']))
            rcov = sum(np.array(covscores) > (special[accession][k] * statistics['coverage']))
            beta = scipy.stats.beta(a+r, b+n-r)
            betacov = scipy.stats.beta(a+rcov, b+n-rcov)

            statistics['beta-p-value-4'] = beta.ppf(0.9999)
            statistics['beta-p-value-3'] = beta.ppf(0.999)
            statistics['beta-p-value-2'] = beta.ppf(0.99)
            statistics['betacov-p-value-4'] = betacov.ppf(0.9999)
            statistics['betacov-p-value-3'] = betacov.ppf(0.999)
            statistics['betacov-p-value-2'] = betacov.ppf(0.99)

            if preview:
                #histcounts, histbins = np.histogram(allbitscores, bins=NUMBINS) 
                #histcounts, histbins = np.histogram(bitscores, bins=NUMBINS) 
                #histcounts, histbins = np.histogram(bitscores, bins=np.arange(np.floor(np.min(bitscores)), np.ceil(np.max(bitscores)+1), 1))
                histcounts, histbins = np.histogram(scores, bins=np.arange(np.floor(np.min(scores)), np.ceil(np.max(scores)+1), 1))
                curarea = np.sum(histcounts * (histbins[1] - histbins[0]))
                normcounts = histcounts / curarea
                histX = (histbins[:-1] + histbins[1:]) / 2
                histY = normcounts
                histYkde = kde.evaluate(histX)

                #X = np.arange(np.min(allbitscores), np.max(allbitscores), (np.max(allbitscores) - np.min(allbitscores)) * 0.001)
                #X = np.arange(np.min(allbitscores), np.max(bitscores), (np.max(allbitscores) - np.min(allbitscores)) * 0.001)
                X = np.arange(np.min(allscores), np.max(scores), (np.max(allscores) - np.min(allscores)) * 0.001)
                Y = kde.evaluate(X)
                plotdata['KDE'] = [X, Y]
                plotdata['Histogram centers'] = [histX, histY]
                plotdata['GEV c unshuf'] = [X, scipy.stats.genextreme.pdf(X, *rawgevparams)]
                plotdata['GEV s unshuf'] = [X, scipy.stats.genextreme.pdf(X, *rawshortgevparams)]
                plotdata['Unshuffled'] = [[special[accession][k]], [scipy.stats.genextreme.pdf(special[accession][k], *rawshortgevparams)]]


                #histcounts, histbins = np.histogram(covbitscores, bins=np.arange(np.floor(np.min(covbitscores)), np.ceil(np.max(covbitscores)+1), 1))
                histcounts, histbins = np.histogram(covscores, bins=np.arange(np.floor(np.min(covscores)), np.ceil(np.max(covscores)+1), 1))
                curarea = np.sum(histcounts * (histbins[1] - histbins[0]))
                normcounts = histcounts / curarea
                histX = (histbins[:-1] + histbins[1:]) / 2
                histY = normcounts
                histYkde = covkde.evaluate(histX)

                #X = np.arange(np.min(covbitscores), np.max(covbitscores), (np.max(covbitscores) - np.min(covbitscores)) * 0.001)
                X = np.arange(np.min(covscores), np.max(covscores), (np.max(covscores) - np.min(covscores)) * 0.001)
                Y = covkde.evaluate(X)
                plotdata['cov KDE'] = [X, Y]
                plotdata['cov Histogram centers'] = [histX, histY]
                plotdata['cov GEV s unshuf'] = [X, scipy.stats.genextreme.pdf(X, *rawshortcovgevparams)]
                #plotdata['cov Unshuffled'] = [[special[accession][k] * statistics['coverage']], [scipy.stats.genextreme.pdf(special[accession][k] * statistics['coverage'], *rawshortcovgevparams)]]
                plotdata['cov Unshuffled'] = [[special[accession][k] * statistics['coverage']], [covkde.evaluate(special[accession][k] * statistics['coverage'])]]

                #plt.plot(X, Y, label='KDE')
                #plt.plot(histX, histY, label='Histogram centers')
                #plt.plot(X, scipy.stats.genextreme.pdf(X, *rawgevparams), label='GEV c unshuf')
                #plt.plot(X, scipy.stats.genextreme.pdf(X, *rawshortgevparams), label='GEV s unshuf')
                #plt.legend()
                #plt.scatter([special[accession][k]], [kde.evaluate(special[accession][k])])
                #plt.show()


        elif 'expect' in k:
            topexpect = 10 ** np.mean(np.log10(bestexpect))
            if special[accession][k] == 0: 
                statistics['expectratio'] = -1
                statistics['medianexpectratio'] = -1
                statistics['log10expectratio'] = -1
                statistics['medianexpect'] = -1
            else: 
                statistics['expectratio'] = topexpect / special[accession][k] 
                #statistics['log10expectratio'] = np.log10(statistics['expectratio'])
                #statistics['medianexpectratio'] = allexpect[len(allexpect)//2] / special[accession][k] 
                #statistics['medianexpect'] = allexpect[len(allexpect)//2]
            statistics['numexpect'] = len(bestexpect)

            #if special[accession][k] > (0.4 * bestexpect[0]):
            #    print(statistics['expectratio'], special[accession][k], 10 ** np.mean(np.log10(bestexpect)), bestexpect, file=sys.stderr)

    statistics['accession'] = accession
    return {'runtime':runtime, 'alignment':special[accession], 'statistics':statistics, 'sequence':sequence}, plotdata

def get_compact_formatted(statsdict):
    out = ''
    #for k in statsdict: print(k, statsdict[k])
    #return out
    try: out += '{queryname}\t{subjectname}\t{identity:0.3f}\t{length}\t\t\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}\t{rank}/{count}\t{swzscore:0.2f}\t{rawpvalue:0.2e}\t{shortpvalue:0.2e}'.format(
        queryname=statsdict['sequence']['_order'][0],
        subjectname=statsdict['sequence']['_order'][1],
        identity=100 * statsdict['alignment']['sw_ident'],
        length=statsdict['alignment']['sw_overlap'],
        qstart=statsdict['sequence']['starts'][0],
        qend=statsdict['sequence']['stops'][0],
        sstart=statsdict['sequence']['starts'][1],
        send=statsdict['sequence']['stops'][1],
        evalue=statsdict['alignment']['sw_expect'],
        bitscore=statsdict['alignment']['sw_bits'],
        rank=statsdict['statistics']['rank'],
        count=statsdict['statistics']['count'],
        swzscore=statsdict['statistics']['swzscore'],
        rawpvalue=statsdict['statistics']['gev-p-value'],
        shortpvalue=statsdict['statistics']['shortgev-p-value'],
    )
    except KeyError: out += '#No alignment found'
    return out


def get_formatted(statsdict):
    out = ''
    out += '########################################\n'
    out += '# Program: {pg_name}\n'.format(**statsdict['runtime'])
    out += '# Rundate: {t} on {host}\n'.format(t=time.strftime("%Y-%m-%dT%H:%M:%SZ%z"), host=socket.gethostname())
    out += '# Commandline: {pg_argv}\n'.format(**statsdict['runtime'])
    out += '# Report_file: stdout\n'
    out += '########################################\n'
    out += '\n'
    out += '#=======================================\n'
    out += '#\n'
    out += '# Aligned_sequences: 2\n'
    out += '# 1: {query}\n'.format(**statsdict['mpsat'])
    out += '# 2: {subject}\n'.format(**statsdict['mpsat'])
    out += '# Matrix: {mp_Parameters}\n'.format(**statsdict['runtime'])
    out += '# Gap_penalty: {pg_open-ext}\n'.format(**statsdict['runtime'])
    out += '#\n'

    if 'sw_overlap' not in statsdict['alignment']: 
        out += '# No alignment found\n'
        return out

    out += '# Length: {sw_overlap}\n'.format(**statsdict['alignment'])

    nident = int(statsdict['alignment']['sw_ident'] * statsdict['alignment']['sw_overlap'])
    pident = '{:0.1%}'.format(statsdict['alignment']['sw_ident']).rjust(4)
    pident = ('(' + pident + ')').ljust(8)
    identval = '{}/{} {}'.format(nident, statsdict['alignment']['sw_overlap'], pident)
    out += '# Identity:{}\n'.format(identval.rjust(20))

    nsim = int(statsdict['alignment']['sw_sim'] * statsdict['alignment']['sw_overlap'])
    psim = '{:0.1%}'.format(statsdict['alignment']['sw_sim']).rjust(4)
    psim = ('(' + psim + ')').ljust(8)
    simval = '{}/{} {}'.format(nsim, statsdict['alignment']['sw_overlap'], psim)
    out += '# Similarity:{}\n'.format(simval.rjust(18))

    #ngaps = int(statsdict['alignment']['sw_
    for k in statsdict['alignment']:
        if 'opt' in k: out += '# Score: {:0.1f}\n'.format(statsdict['alignment'][k])
    #for k in statsdict['alignment']:
    #    if 'expect' in k: out += '# Expect: {:0.1e}\n'.format(statsdict['alignment'][k])
    for k in statsdict['statistics']:
        if 'coverage' in k: out += '# Coverage: {:0.1%}\n'.format(statsdict['statistics'][k])
    
    
    out += '#\n'
    out += '#\n'
    out += '#=======================================\n'
    out += '\n'

    qname = statsdict['sequence']['_order'][0]
    sname = statsdict['sequence']['_order'][1]
    qstr = statsdict['sequence'][qname]
    sstr = statsdict['sequence'][sname]
    cons = statsdict['sequence'][statsdict['sequence']['_order'][2]].replace(':', '|').replace('.', ':').replace('-', ' ')
    qresi = 1
    sresi = 1
    length = 50
    if length:
        for start in range(0, len(qstr), length):
            qseg = qstr[start:start+length]
            sseg = sstr[start:start+length]

            if len(qseg) < len(sseg): qseg += '-' * (len(sseg) - len(qseg))
            elif len(sseg) < len(qseg): sseg += '-' * (len(qseg) - len(sseg))
            nextqresi = qresi + len(qseg) - qseg.count('-')
            nextsresi = sresi + len(sseg) - sseg.count('-')

            endqresi = qresi if qseg.count('-') == length else nextqresi - 1
            endsresi = sresi if sseg.count('-') == length else nextsresi - 1

            out += '{}{} {}{}\n'.format(qname.ljust(13)[:13], str(qresi).rjust(7), qseg, str(endqresi).rjust(7))
            out += '{}{}\n'.format(''.ljust(21)[:21], cons[start:start+length])
            out += '{}{} {}{}\n'.format(sname.ljust(13)[:13], str(sresi).rjust(7), sseg, str(endsresi).rjust(7))
            qresi = nextqresi
            sresi = nextsresi
            out += '\n'
    else:
        out += '{}{}\n'.format(qname.ljust(21)[:21], qstr)
        out += '{}{}\n'.format(''.ljust(21)[:21], cons)
        out += '{}{}\n'.format(sname.ljust(21)[:21], sstr)

    out += '#\n'
    out += '# Rank: {rank}/{count}\n'.format(**statsdict['statistics'])
    out += '# Padding: {npad}\n'.format(**statsdict['statistics'])
    out += '# S-W z-score: {swzscore:0.1f} ({swmean:0.1f} +/- {swstd:0.1f})\n'.format(**statsdict['statistics'])
    #out += '# Bits z-score: {bitszscore:0.1f} ({bitsmean:0.1f} +/- {bitsstd:0.1f})\n'.format(**statsdict['statistics'])
    if statsdict['statistics'] == -1: 
        pass
        #out += '# Expect ratio: Too high (identical sequences?)\n'
        #out += '# Log10-expect ratio: Too high (identical sequences?)\n'
        #out += '# Median expect ratio: Too high (identical sequences?)\n'
    else: 
        #out += '# Expect ratio (top {numexpect}): {expectratio:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Log10-expect ratio (top {numexpect}): {log10expectratio:0.1f}\n'.format(**statsdict['statistics'])
        #out += '# Median expect ratio: {medianexpectratio:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Median expect: {medianexpect:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Gaussian KDE p-value (no adjustment): {p-value:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Histogram-KDE mean absolute error: {histkde_mae:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Histogram-KDE mean squared error: {histkde_mse:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Gaussian KDE p-value [unshuf, 2*maxval]: {p-value-2x:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Raw GEV p-value with unshuf [unshuf, inf]: {gev-p-value:0.2e}\n'.format(**statsdict['statistics'])
        out += '# Gaussian KDE p-value (coverage adjusted): {shortcov-p-value:0.2e}\n'.format(**statsdict['statistics'])
        out += '# Raw GEV p-value (no adjustment): {shortgev-p-value:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Raw GEV p-value (coverage adjusted): {shortcovgev-p-value:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Inferred P-value (0.9999): {beta-p-value-4:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Beta p-value (0.999 conf, no adjustment): {beta-p-value-3:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Inferred P-value (0.99): {beta-p-value-2:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Inferred P-value with coverage adjustment (0.9999): {betacov-p-value-4:0.2e}\n'.format(**statsdict['statistics'])
        out += '# Beta p-value (0.999 conf, coverage adjusted): {betacov-p-value-3:0.2e}\n'.format(**statsdict['statistics'])
        #out += '# Inferred P-value with coverage adjustment (0.99): {betacov-p-value-2:0.2e}\n'.format(**statsdict['statistics'])

    return out

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--flags', help='Flags to pass on to ssearch as a quoted string. -E, -m, and -k will be ignored.')
    parser.add_argument('-c', '--count', type=int, default=1000, help='Number of shuffles to perform (default:1000)')
    parser.add_argument('-e', type=int, default=10, help='Number of top shuffled alignments to compare the unshuffled alignment against (default:10)')
    parser.add_argument('-i', action='store_true', help='Prompt with filenames to allow manual inspection')
    #parser.add_argument('-m', '--mode', default='swzscore', help='Statistic to use (\033[1mswzscore\033[0m, bitzscore, rank)')
    parser.add_argument('-o', '--outdir', help='Where to write intermediate files (default: nowhere)')

    parser.add_argument('-q', '--quiet', action='store_true', help='Don\'t print anything to stdout')

    #mergeparams = parser.add_mutually_exclusive_group()
    mergeopts = parser.add_argument_group("Mergedict tuning arguments. Mutually exclusive with --stms and --sub-states")
    mergeopts.add_argument('-s', action='store_true', help='Split loop residues by in/out, i.e. disable I-O and i-o merging')
    mergeopts.add_argument('-t', action='store_true', help='Split residues by loop/tail, i.e. disable I-i and O-o merging')

    parser.add_argument('--stms', nargs='*', help='Space-separated ranges to use as TMSs instead of HMMTOP\'s predictions for the subject protein. Mutually exclusive with -s/-t and --sub-states')
    parser.add_argument('--sub-states', help='Use this state string for the subject. Bypasses --stms and mergedict-related arguments. Mutually exclusive with -s/-t and --stms')

    parser.add_argument('--binsize', type=float, default=None, help='Bin size for histogram (default: auto)')
    parser.add_argument('--bandwidth', type=float, default=None, help='Bandwidth for (default: auto (Scott\'s method)')

    parser.add_argument('--show', action='store_true', help='Plot the KDE and the histogram')
    parser.add_argument('--resume', default=False, nargs='?', help='Run statistics on an existing LSAT output directory or SSEARCH results generated in the machine-readable format (-m 10). If used without an argument, --resume will attempt to use the value of --outdir')
    parser.add_argument('--compact', action='store_true', help='Enable compact output for large-scale analyses')
    parser.add_argument('--gsat', action='store_true', help='Run in GSAT mode by ignoring all TMSs. Takes precedence over --stms')
    parser.add_argument('--pad-unalignable', action='store_true', help='Copy and use the worst score for all unalignable sequences')

    parser.add_argument('query', nargs='?')
    parser.add_argument('subject', nargs='?')

    figopts = parser.add_argument_group('Figure tuning arguments.')
    figopts.add_argument('--save', nargs='?', default=-1, help='Save figure. Defaults to OUTDIR/plot.png if given without arguments')
    figopts.add_argument('--dpi', type=int, default=None, help='Figure DPI')
    figopts.add_argument('--height', type=int, default=None, help='Figure height')
    figopts.add_argument('--width', type=int, default=None, help='Figure width')

    args = parser.parse_args()

    #if args.mode not in ('swzscore', 'bitzscore', 'rank'): error('Invalid mode "{}"'.format(args.mode))
    if args.resume is False and args.subject is None:
        parser.print_usage()
        exit(1)
    if sum([args.stms is not None, args.sub_states is not None, (args.s or args.t)]) >= 2:
        parser.error('argument --stms, --sub-states, and -s/-t are mutually exclusive')
        exit(1)

    main(args.query, args.subject, count=args.count, outdir=args.outdir, flags=args.flags, split_inout=args.s, split_looptail=args.t, interactive=args.i, expectsmooth=args.e, quiet=args.quiet, binsize=args.binsize, bandwidth=args.bandwidth, preview=args.show, resume=args.resume, compact=args.compact, stms=args.stms, gsat=args.gsat, pad_unalignable=args.pad_unalignable, dpi=args.dpi, height=args.height, width=args.width, save=args.save)
