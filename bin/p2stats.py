#!/usr/bin/env python
import protocol2
import sys, os
import numpy

usage = 'p2stats.py <min_SD> <overlap> <tms,tms>'
if len(sys.argv) is not 4:
    print usage
    quit()

minsd,overlap,restrict = sys.argv[1:4]
protocol = protocol2.Compare()
if os.path.exists('.'+protocol.tssout) is False or  os.path.exists('.'+protocol.gsatout) is False:
    print 'CD to an existing P2 outputfolder'
    quit()
protocol.outdir = '.'
protocol.tssearch()
protocol.run_global_alignments()
sr,tr = restrict.split(',')
for result in protocol.globaldata:
    subject,target,score = result['subject_annotation'],result['target_annotation'],result['gsat_score']
    try:
        if score >= int(minsd):
            if bool(int(sr)) or bool(int(tr)):
                if result['subject_tmc'] != int(sr) or result['target_tmc'] != int(tr):
                    continue
            protocol.calculate_overlap(subject,target,score,float(overlap))
    except:
        pass

for pairs, data in protocol.stats.items():
    print pairs,
    print numpy.average(data),
    print numpy.std(data)
