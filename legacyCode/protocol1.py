#!/usr/bin/env python

'''
   A more formal solution to remote PSI BLASTS. NCBI
   Has changed their format, and broke the PHP version.
   This Mecahnized edition should be hold up against
   Most minor changes in the future.

   Designed by Vamsee Reddy. Part of the BioV Suite.
'''

import mechanize
import sys,re
import cookielib
import tempfile
from time import sleep
from Bio.Blast import NCBIXML
from Bio import SeqIO

class psi:

    def __init__(self):

        # PSI Blast Settings
        self.expect = 10
        self.query = None
        self.results = 500
        self.minaln = 60

        # Internal Browser settings
        self.gis = []
        self.br = mechanize.Browser(factory=mechanize.RobustFactory())
        cj = cookielib.LWPCookieJar()
        self.br.set_cookiejar(cj)
        self.br.set_handle_equiv(True)
        self.br.set_handle_redirect(True)
        self.br.set_handle_referer(True)
        self.br.set_handle_robots(True)
        #self.br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
        self.br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) \
        Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

    def init_blast(self):
        self.br.open('http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=blastp')
        self.br.select_form(nr=0)
        self.br.form.set_all_readonly(False)
        c = self.br.response().read()
        self.br.form['SELECTED_PROG_TYPE'] = 'psiBlast'
        self.br.form['RUN_PSIBLAST'] = 'on'
        self.br.form['BLAST_PROGRAMS'] = ['psiBlast']
        self.br.form['QUERY'] = self.query
        self.br.form['I_THRESH'] = str(self.expect)
        self.br.form['EXPECT'] = str(self.expect)
        self.br.form['WORD_SIZE'] = ['3']
        self.br.form['DESCRIPTIONS'] = str(self.results)
        self.br.form['ALIGNMENTS'] = str(self.results)
        self.br.form['NUM_OVERVIEW'] = str(self.results)
        self.br.form['MAX_NUM_SEQ'] = [str(self.results)]
        self.br.submit()
        self.loop()

    def wait(self):
        text = re.compile('This page will be automatically updated in <b>(\d+)<\/b> seconds?')
        time = text.search(self.br.response().read())
        if bool(time):
            seconds = int(time.groups()[0])
            secs = seconds if seconds < 15 else 10
            print 'Waiting %i Seconds' %secs
            sleep(secs)
            return True
        return False

    def loop(self):
        if self.wait():
            self.br.select_form(nr=0)
            self.br.submit()
            self.loop()
        return

    def fetch_results(self):
        for link in self.br.links(url_regex='FORMAT_TYPE=XML'):
            req = self.br.click_link(link)
            ### inserted here from below
            (path,obj) = self.br.retrieve(req)
            gis = []
            xml = open(path,'r')
            #print xml.read()
            #xml.seek(0)
            for blast in NCBIXML.parse(xml):
                for desc in blast.descriptions:
                    gi = re.search(r'^gi\|(\d+)\|',str(desc)).groups()[0]
                    longEnough = 0
                    for aln in blast.alignments:
                        for hsp in aln.hsps:
                            testLn = hsp.query.replace('-','')
                            if len(testLn) >= self.minaln:
                                longEnough += 1
                    if longEnough > 0:
                        gis.append(gi)
            #self.gis = list(set(gis))
            if len(gis) > 0:
                if len(self.gis) > 0:
                    gis.extend(self.gis)
                    self.gis = list(set(gis))
                else:
                    self.gis = list(set(gis))
            ### end of insertion from below
            break
        # moved above
        #(path,obj) = self.br.retrieve(req)
        #gis = []
        #xml = open(path,'r')
        #for blast in NCBIXML.parse(xml):
        #    for desc in blast.descriptions:
        #        gi = re.search(r'^gi\|(\d+)\|',str(desc)).groups()[0]
        #        gis.append(gi)
        #self.gis = list(set(gis))

    def iterate(self):
        #results = self.br.response().read()
        #results = open('out.html','r').read()
        #results = results.replace('<! --','<!--')
        #myres = tempfile.NamedTemporaryFile(suffix='.html', delete=False)
        #print myres.name
        #myres.write(results)
        #myres.seek(0), myres.flush()
        #self.br.open_local_file(myres.name)
        self.br.select_form(nr=3)
        #self.br.form.action='http://blast.ncbi.nlm.nih.gov/Blast.cgi'
        self.br.submit()
        self.loop()


###### here the program runs
if __name__=='__main__':

    import blast
    import optparse
    import shutil,os
    import subprocess
    thisprogname = os.path.basename(__file__)
    desc = '''
    Welcome to {0}.
    This tool will run a PSI Blast with iterations, collect results,
    remove redundant/similar sequences annotate, tabulate, & count TMSs.
    Developed by Vamsee Reddy :: Part of the BioV Suite.
    '''.format(thisprogname)
    desc = " ".join(desc.split())
    opts = optparse.OptionParser(description=desc,version=2.0)
    opts.add_option('-q',
                    action='store',
                    dest='query',
                    type='string',
                    default=None,
                    help='Gi/Accession/Sequence to BLAST.')
    opts.add_option('-i',
                    action='store',
                    dest='iterate',
                    type='int',
                    default=1,
                    help='Number of additional iterations to perform. (1)')
    opts.add_option('-n',
                    action='store',
                    dest='number',
                    type='int',
                    default=500,
                    help='Number of results to fetch each round (500)')
    opts.add_option('-e',
                    action='store',
                    dest='expect',
                    type='float',
                    default=0.005,
                    help='E-Value cutoff (0.005)')
    opts.add_option('-c',
                    action='store',
                    dest='cutoff',
                    type='float',
                    default=0.8,
                    help='CD-HIT threshold. From 0.4 - 1 (0.8)')
    opts.add_option('-o',
                    action='store',
                    dest='outdir',
                    type='string',
                    default='p1out',
                    help='Output folder (p1out)')
    opts.add_option('--tms',
                    action='store_true',
                    dest='tms',
                    default=False,
                    help='Include this flag to tabulate TMS stats.')
    opts.add_option('--min',
                    action='store',
                    dest='minaln',
                    default=60,
                    type='int',
                    help='Minimum sequence length to retrieve')
    opts.add_option('--max',
                    action='store',
                    dest='maxaln',
                    default=0,
                    type='int',
                    help='Maximum sequence length to retrieve')
    (cli,args) = opts.parse_args()
    protocol   = psi()
    increments = [10,50,100,250,500,1000,5000,10000,20000]
    if cli.number not in increments:
        print "NCBI can only retrieve sequences in these increments:\n"+', '.join([str(i) for i in increments])
        quit()
    protocol.results = cli.number
    protocol.expect  = cli.expect
    protocol.minaln  = cli.minaln
    if cli.query is None:
        opts.print_help()
        print '\n--------\n# User Input Options'
        # Simple approach for students :
        protocol.query = raw_input('Accession/Gi/Seq? ')
        hmmtop_bool = raw_input('Count TMSs? (y/N) ')
        outdir = raw_input('Output path? ')
        protocol.init_blast()
        protocol.fetch_results()
        print 'Found %i results' %len(protocol.gis)
        if bool(cli.iterate):
            for i in range(cli.iterate):
                print 'Running iteration round #%i'%(i+1)
                protocol.iterate()
                protocol.fetch_results()
                newFound = len(protocol.gis) - lastFound
                print 'Found %i new results'%newFound
                ### added break if no more results are found
                if newFound < 1:
                    print '  no more iterations possible'
                    break
                else:
                    lastFound += newFound
        cutoff = raw_input('CD-HIT threshold (0.4 - 1)? ')
    else:
        outdir,cutoff,hmmtop_bool = cli.outdir,cli.cutoff,''
        protocol.query = str(cli.query)
        protocol.init_blast()
        protocol.fetch_results()
        print 'Found %i results'%len(protocol.gis)
        lastFound = len(protocol.gis)
        if bool(cli.iterate):
            for i in range(cli.iterate):
                print 'Running iteration round #%i'%(i+1)
                protocol.iterate()
                protocol.fetch_results()
                newFound = len(protocol.gis) - lastFound
                print 'Found %i new results'%newFound
                ### added break if no more results are found
                if newFound < 1:
                    print '  no more iterations possible'
                    break
                else:
                    lastFound += newFound
    ### there's a problem here if iterations erase gis list:
    blast = blast.tools()
    blast.gis = protocol.gis
    # blast.build_xml()

    if os.path.exists(outdir) is False:
        os.mkdir(outdir)
    blast.build_raw_fasta()
    shutil.copy(blast.raw_fasta.name,outdir+'/results')

   ### shutil.copy(blast.xml_file.name,outdir+'/results')
    cutoff = float(cutoff) if cutoff != '' else int(1)
    if bool(cutoff):
        ### cmd = 'make_table5.pl -i "%s" -s 10 -l 999999 -c %f' %(outdir+'/results',cutoff)
        if cutoff > 0.7:
            cmd = 'cd-hit -i %s -o %s -d 50 -c %f'%(outdir+'/results',outdir+'/results.faa',cutoff)
        if cutoff > 0.6 and cutoff <=0.7:
            cmd = 'cd-hit -i %s -o %s -n 4 -c %f'%(outdir+'/results',outdir+'/results.faa',cutoff)
        if cutoff > 0.5 and cutoff <= 0.6:
            cmd = 'cd-hit -i %s -o %s -n 3 -c %f'%(outdir+'/results',outdir+'/results.faa',cutoff)
        if cutoff < 0.5 and cutoff >= 0.4:
            cmd = 'cd-hit -i %s -o %s -n 2 -c %f'%(outdir+'/results',outdir+'/results.faa',cutoff)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout,stderr = process.communicate()


    if hmmtop_bool.lower() == 'y' or hmmtop_bool.lower() == 'yes' or cli.tms is True:
        import hmmtop,tmstats
        (stats,tms) = tmstats.calculate(outdir+'/results.faa','Protocol1 TMS Stats',outdir+'/tms_stats.png')
        tab_file = open(outdir+'/results.tab','r')
        lines = tab_file.read().strip().split('\n')
        tab_file.close()
        tab_file = open(outdir+'/results.tab','wb')
        for i,line in enumerate(lines):
            line = line.split('\t')
            tmc = len(tms[line[0]]) if line[0] in tms else 0
            line[1:1]=[str(tmc)] if i > 0 else ['TMSs']
            tab_file.write('\t'.join(line)+'\n')
        tab_file.close()

    # Removes specified long/short seqs. Added 12/13/11 -Vamsee
    if bool(cli.minaln) or bool(cli.maxaln):
        valids = []
        mymax = cli.maxaln if bool(cli.maxaln) else 99999999
        fastas = SeqIO.parse(open(outdir+'/results.faa'),'fasta')
        [valids.append(i) for i in list(fastas) if len(i.seq) >= cli.minaln and len(i.seq) <= mymax]
        handle = open(outdir+'/results.faa','wb')
        for i in valids:
            SeqIO.write(i,handle,'fasta')

