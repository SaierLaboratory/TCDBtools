#!/usr/local/anaconda3/bin/python3


'''
Hmmtop module will allow you to import FASTA sequences, HMMTOP all of
them, and retrieve TMS Ranges for any sequence by symbol. Results are
preserved throughout a session, so you can Pass it to other programs
instead of running it over and over again. Supports multiple
Libraries, so you can load several FASTA files that share the same
symbol safely.

ht = hmmtop.tools() >> ht.add_library(name,path) >> \
   ht.scan_libraries() (runs hmmtop).

get results from self.results in a very cleary dictionary format.
'''
import sys
import tempfile
import subprocess
import hashlib
import re,settings,os,pickle
from ProjectBio import ParseDefline

class tools:

    def __init__(self):
        self.results = {}
        self.libraries = {}
        self.pop=[]

    def add_library(self,name,path):
        if name in self.libraries:
            raise IOError("This library already exists")
            return
        self.libraries[name] = path
        return

    def scan_file(self,fasta_file):
        hmtool = 'hmmtop' if 'hmmtop' not in os.environ else os.environ['hmmtop']
        cmd = "%s -if=%s -sf=FAS -is=pseudo -pi=spred" %(hmtool,fasta_file)

        handle = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)

        results = handle.communicate()[0]
        tms = re.compile("\s+(IN|OUT)\s+([0-9]+)\s+([^\n]+)?")
        ranges = re.compile("(\d+)\s+(\d+)")
        name = re.compile("^>HP:\s+\d+\s+(\S+)")
        
        lines = [i for i in results.rstrip().split('\n') if len(i) > 1]
        (res,symbols)=[],[]

        for pos,i in enumerate(lines):

            #Ignore matches with 0 TMSs
            if re.search('(IN|OUT)\s+0', i):
                continue

            line = tms.search(i).groups()

            keys = {}
            for x in enumerate(ranges.finditer(line[2])):
                keys.setdefault(x[0]+1,[]).append(int(x[1].groups()[0]))
                keys.setdefault(x[0]+1,[]).append(int(x[1].groups()[1]))
            #symbol used to be name.search(i).groups()[0].replace('>','')
            symbols.append(ParseDefline(name.search(i).groups()[0]).id)
            res.append(keys)
        return res, symbols

    def scan_libraries(self):
        if ( 'BIOV_DEBUG' in os.environ and
             os.environ['BIOV_DEBUG'].lower() == 'true' ):
            hashtable = []
            for name,fn in self.libraries.items():
                hashtable.append(open(fn,'rb').read())
            debugfile = '/tmp/'+hashlib.md5(''.join(hashtable)).hexdigest()+'.hmmtop'
            if os.path.exists(debugfile):
                self.results = pickle.load(open(debugfile,'r'))
                return
        
        #Python 3 does not support pop operations while interating inside a 
        #dictionary anymore. So, a copy should be made with:
        #              "list(self.libraries.items())"
        for name,file in list(self.libraries.items()):
            (result,symbols) = self.scan_file(file)
            
            keys = dict(zip(symbols,result))
            self.results[name] = keys
            self.libraries.pop(name)

        if ( 'BIOV_DEBUG' in os.environ and
             os.environ['BIOV_DEBUG'].lower() == 'true' ):
            pickle.dump(self.results,open(debugfile,'wb'))

        return
