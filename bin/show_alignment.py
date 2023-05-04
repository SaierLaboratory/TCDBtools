#!/usr/bin/env python
from Bio import SeqIO
import os,sys,gsat,re,tempfile,hmmtop,copy
get = raw_input("Input: ")
oget=get
maxflank = 20
try:
    subjects = open("subject.faa",'r')
    targets = open('target.faa','r')

    get = get.strip().split("\t")
    gs = get[1]
    gt = get[3]    

except:
    data = re.compile('(\w+)\t([0-9,-]+) && ([0-9,-]+)')
    dat = data.findall(get)[0]
    gs = "%s TMS #(%s)"%(dat[0],dat[1].replace('-',','))
    gt = "%s TMS #(%s)"%(dat[0],dat[2].replace('-',','))
    subjects = open("./"+dat[1]+'/subject.faa','r')
    targets = open("./"+dat[1]+'/target.faa','r')

finally:

    subjects = SeqIO.parse(subjects,'fasta')
    targets = SeqIO.parse(targets,'fasta')

    sub = {}
    tar = {}

    [sub.setdefault(i.description,str(i.seq)) for i in subjects]
    [tar.setdefault(x.description,str(x.seq)) for x in targets]
    s = sub[gs]
    t = tar[gt]
    print s,t
    gsat.compare(s,t,'global',100,8,2,True)

    if raw_input("Optimize Alignment? (y/n) :").lower() != "y":
        exit()

    consec = re.compile('([0-9,-]+) && ([0-9,-]+)')
    consecv = re.compile('\(([0-9,-]+)\).+\(([0-9,-]+)\)')
    tmss = consec.search(str(oget)).groups() if bool(consec.search(str(oget))) is True else consecv.search(str(oget)).groups()

    ff = raw_input('Where is your original FASTA File? :').strip()
    fasta = open(ff,'r')
    fasta = SeqIO.parse(fasta,'fasta')
    fasta = SeqIO.to_dict(fasta)

    asequence = fasta[gs.split(' ')[0]]
    bsequence = fasta[gt.split(' ')[0]]

    optfile = tempfile.NamedTemporaryFile()
    SeqIO.write(asequence,optfile,'fasta')
    SeqIO.write(bsequence,optfile,'fasta')
    optfile.seek(0)
    ht = hmmtop.tools()
    ht.add_library('OPT',optfile.name)
    ht.scan_libraries()
    TMSR = ht.results['OPT']

    def expand(seq,tms,num,select):
        in_use = []
        tmss = copy.deepcopy(tms)
        [in_use.extend(range(i[0],i[1]+1)) for i in tmss.values()]
        maxlen = len(seq)
        loops = []
        for tms,residues in tmss.items():
            for extend in range(num):
                left = residues[0] - 1
                right = residues[1] + 1
                if left >= 1 and left not in in_use:
                    in_use.append(left)
                    residues[0] = left
                if right <= maxlen and right not in in_use:
                    in_use.append(right)
                    residues[1] = right
            loops.append(seq[residues[0]-1:residues[1]])
        res = ''.join([loops[x-1] for x in select])
        return res

    gsats = []
    aselect = [int(i) for i in re.split(",|-",tmss[0])]
    bselect = [int(i) for i in re.split(",|-",tmss[1])]

    for xx in range(maxflank+1):
        print "-",
        a = expand(str(asequence.seq),TMSR[asequence.id],xx,aselect)
        b = expand(str(bsequence.seq),TMSR[bsequence.id],xx,bselect)
        gs = gsat.compare(a,b,'global',500,8,2,False)
        print xx, gs.zscorep
        gsats.append(gs)
    gsats.sort(key=lambda x:x.zscorep,reverse=True)
    gsats[0].outfile.seek(0)
    print gsats[0].outfile.read()
    print "\n----------------------------\n OPTIMAL Z-SCORE IS : %i"%gsats[0].zscore
    print "%s :: %s \n\n%s :: %s"%(asequence.id,gsats[0].AA,bsequence.id,gsats[0].BB)



