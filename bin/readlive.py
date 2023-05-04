#!/usr/bin/env python
import sys,os,re
consec = re.compile('([0-9,-]+) && ([0-9,-]+)')
consecv = re.compile('\(([0-9,-]+)\).+\(([0-9,-]+)\)')
from time import sleep
while 1:
    handle = open(sys.argv[1],'r')
    data = []
    handle.seek(0)
    for i in handle.readlines():
        line = i.strip().split("\t")
        data.append(line)
    data.sort(key=lambda x:(100-float(x[-2]),float(x[-1])))
    for v in data:
        thisline = "\t".join(v)
        score = float(v[-1])
        if score < int(sys.argv[3]):
            continue
        tmss = consec.search(thisline).groups() if bool(consec.search(thisline)) is True else consecv.search(thisline).groups()
        if sys.argv[2] == 'consec':
            start = int(re.split(',|-',tmss[0])[-1])
            end = int(re.split(',|-',tmss[1])[0])
            if start + 1 == end:
                print thisline
        if sys.argv[2] != 'all' and sys.argv[2] !='consec':
            if tmss[1] ==sys.argv[2]:
                print thisline
        if sys.argv[2] == 'all':
            print thisline

    print "#%i Results... "%len(data)
    what = raw_input("Push Enter to refresh...")
    os.system('clear')

