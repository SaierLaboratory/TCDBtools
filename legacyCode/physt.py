#!/usr/bin/env python
'''
  phyST
  This program can take one or more TCDB family codes (i.e. 1.A.1) as input and runs PSI-BLAST for all proteins in each family.
  The taxanomic classification for each BLAST hit is extractcted and the protein composition of each detected phyla within the query family(-ies) is reported.
  PhyST can generate the FASTA sequences of top scoring hits, which can be used to run other programs (e.g. multiple alignment, hydropathy analyses, etc).
'''
# Written by Hari Krishnan, Larry Chau
# lchau@ucsd.edu - larrymchau@gmail.com
# hkkrishn563@gmail.com - hkkrishn@ucsd.edu

import ntpath
import mechanize
import sys,re,os
import urllib2
import tempfile
from Bio import Entrez
from bs4 import BeautifulSoup
from time import sleep
from urllib2 import Request,urlopen,URLError, HTTPError
import httplib
import cookielib
import math
import multiprocessing as mp
import argparse
import ctypes
import pprint
import numpy as np
import math
import pdb
import time

#Specify the user's email when using Entrez from Bio package
Entrez.email = "hkkrishn@ucsd.edu"
    
#Globals, working data
interval = 0
br = mechanize.Browser()
tfiles = []
lo_cutoff = 50
eValue = '0.0001'
alignment = ''
word_size = '3'
start_time = 0
stop_after = 800 #seconds
timed_out = False
tmstrue = True
datatrue = True

def browser_init():
    global br
    cj = cookielib.LWPCookieJar()
    br.set_cookiejar(cj)
    br.set_handle_equiv(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
    br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) \
    Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

#Code copied from protocol1.py
def wait():
    global interval
    global br
    
    text = re.compile('This page will be automatically updated in <b>(\d+)<\/b> seconds?')
    time = text.search(br.response().read())
    if bool(time):
        # for x in range (0,40):
        seconds = int(time.groups()[0])
        secs = seconds if seconds < 15 else 10
        if(interval >= 15):
           interval = 0
           print '   waited 15 seconds...'
        interval = interval + secs
        sleep(secs)
        return True
    return False

#Code copied from protocol1.py
def loop():
    global br
    global start_time
    global timed_out
    global stop_after
    
    current_time = time.time()
    if(current_time - start_time) >= stop_after:
        print "    There seems to be a connection Timeout."
        timed_out = True
        return
    if wait():
        br.select_form(nr=0)
        br.submit()
        loop()
        return

def find_optimal(proteinCode):
    global br
    global alignment
   
    br.open('http://www.tcdb.org/progs/blast.php')
    br.select_form(nr=1)
    br.form['BLAST'] = ['blastp']
    br.form['SEQUENCE'] = alignment
    br.form['EXPECT'] = ['1000']
    br.form['DESCRIPTIONS'] = ['200']
    response1 = br.submit().read()
    response1 = response1.split('<pre>')[2]
    text = re.split('<a href="/search/result.php\?tc=\d+.\D+.\d+">',response1)
    del text[0]

    for i in range(0,len(text)):
         family = text[i].split('</a>')[0]
        #only match substring to tcdb family code
         if proteinCode != family[:len(proteinCode)]:
            if (i != 0): optimEVal = text[i-1].split('</a>')[2].split('<a href="')[0].strip()
            break

    nums = optimEVal.split('e')
    if(1 < len(nums)):
        if(re.match('\W*',nums[0])): nums[0] = '1'
        optimEVal = float(nums[0])*(10**float(nums[1]))
    else: optimEVal = float(nums[0])
    optimEVal = '%.9f' %optimEVal    
    return optimEVal

def blast(proteinCode,fileFlag):
    global eValue
    global alignment
    global br
    global word_size
    global timed_out

    alignment = ""

    if(fileFlag):

        if os.path.exists(proteinCode):
            f = open(proteinCode)
            all_seqs = f.read()
            split_seqs = all_seqs.split('\n>')
                
            for seq in split_seqs:
                    
                if not bool(re.search('^>', seq)):
                    alignment = ">" + seq
                else:
                    alignment = seq
                
                if alignment:
                    break
                
            f.close()

    else:
        browser_init()

        #Puts the family code into tcdb.org and search
        br.open('http://www.tcdb.org/')
        br.select_form(nr=0)
        br.form.set_all_readonly(False)
        br.form['query'] = proteinCode
        br.submit()

        if (len(proteinCode.split('.')) < 4):
          #Clicks the link containing text "View Proteins beloning to blah"
          link = br.click_link(text_regex = "View Proteins belonging to: ")
          br.open(link)

        #Click on the first subfamily on the subfamily list.
        cnt = 0
        while True:
          link = br.click_link(text_regex = proteinCode, nr= cnt)
          response = br.open(link)
          #  Incase that it is possible to not entering a protein's info page
          #  after clicking "View proteins", skip the first subfamily and
          #  go to the next one, etc.
          try:
              #The expected FASTA page is the link with text "FASTA
              #  formatted sequence", which contains url regex "fasta.php"
              link = br.click_link(url_regex = "fasta.php")
              break
          except mechanize._mechanize.LinkNotFoundError:
              #If the page does not contain "fasta.php", skip to the next
              #  subfamily
              br.back()
              cnt = cnt + 1

        print link

        #click into the FASTA fornatted sequence, then split texts to
        #  extract the string containing only alignment sequence
        sourcePage = br.open(link)
        keyLines = sourcePage.read().split('<PRE>')[1]
        keyLines = keyLines.split('</PRE>')[0]
        keyLines = keyLines.split('\n')
        del keyLines[0]
        for row in keyLines:
          alignment  = alignment + row

        optimEVal = find_optimal(proteinCode)

        print '    Estimate Optimal E-Value (TCDB):',optimEVal,'(e.g. ',float(optimEVal),')'
        print '    Using:',eValue
        eValue = eValue.replace(' ', '')

        #Go to NCBI blast page, enter the alignment found above, select
        #  "psi-blast" and "5000" max results", then blast
        req = Request('http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome')
        try:
            response = urlopen(req)
            print "    Connecting to BLAST server..."
            br.open('http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome')
            br.select_form(nr=0)
            br.form.set_all_readonly(False)
            br.form['SELECTED_PROG_TYPE'] = 'psiBlast'
            br.form['RUN_PSIBLAST'] = 'on'
            br.form['BLAST_PROGRAMS'] = ['psiBlast']
            br.form['QUERY'] = alignment
            br.form['I_THRESH'] = eValue
            br.form['MAX_NUM_SEQ'] = ['10000']
            br.form['WORD_SIZE'] = [word_size]
            br.submit()
            print "    Connected!"
            print "    Blasting Off..."
            loop()
            if( not timed_out ):
                print "    Done.\n"
            return
        except:
            print "    Sorry, The BLAST server is not  responding. Please check your network connectivity."
            timed_out = True
            return

# print "Oops! It seems that the server is not responding I'll try again"

# To explain this code a little more, it is easier to write the response
# to a html file and submit the file using mechanize parameters, for
# there does not seem to be any inherent support to follow a response
# without a link. In this case, the blast returns a custom source page
# through the CGI script. The only way to submit the content again
# is to hard code the response. As we can see done here.
def iterate():
    global br
    global tfiles
    
    results = br.response().read()
    results = results.replace('<! --','<!--')
    
    myres = tempfile.NamedTemporaryFile(mode='w+t',suffix='.html', delete=False)
    myres.write(results)
    myres.seek(0), myres.flush()
    br.open_local_file(myres.name)
    myres.close()
    
    # find and select form
    formcount=0
    for form in br.forms():
        if 'name' in form.attrs:
            if form.attrs['name'] == 'overview0':
                br.form = form
                break
                formcount=formcount+1

    br.select_form(nr=formcount)
    br.form.action='http://blast.ncbi.nlm.nih.gov/Blast.cgi'
    br.submit()
    loop()
    
    tfiles.append(myres.name)
    return

def close():
    global tfiles
    for item in tfiles:
        os.remove(item)

def process_data(process, keyLines, dataQueue, phylumQueue, sequencesQueue, invalidQueue, cutoff):
    
    invalidProteins = 0
    
    minLength = int((len(keyLines)*0.25)*process)
    maxLength = int((len(keyLines)*0.25)*(process+1))
    counter = 0
    
    #For each link containing the accession number as link
    for item in keyLines[minLength:maxLength]:
        #Extract the e-value
        eValue = (item.split('<td>')[4]).split(' </td')[0]
        eValue = float(eValue.strip('\n'))

        #Extract the accession number
        accessionNum = (item.split('Show report for ')[1]).split('"')[0]

        #Use efetch() to fetch info of this accession number, from database
        #  "protein", return type as GenPept flat file, return mode as xml
        #  More info: http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        try:
            handle = Entrez.efetch(db="protein", id=accessionNum, rettype="gp", retmode="xml")
        except urllib2.HTTPError:
            #Incase that this accession number is not in "protein" database
            #  (e.g. 3D structure as "Chain A, Potassium Channel (Kcsa)
            #  Full-Length Fold", accession "1F6G_A"
            continue

        #read the xml, and extract info: sequence length, domain, phylum,
        #  and protein name. Filter the proteins by lower cutoff length
        record = Entrez.read(handle)

        try:
            phylum = record[0]['GBSeq_taxonomy'].split(';')[1]
        except IndexError:
            #Incase that xml does not contain taxonomy info, skip to next
            #  protein. Cause unknown.
            continue

        seqLen = int(record[0]['GBSeq_length'])
        
        query_cover = (item.split('<td>')[3]).split('%')[0]
        query_cover = int(query_cover.strip('\n'))
        if query_cover < cutoff:
            invalidQueue.put(1)
            continue

        try:
            length = len(record[0]['GBSeq_other-seqids'])
            for i in range(0,length):
                if('gi' == record[0]['GBSeq_other-seqids'][i].split('|')[0]):
                    giNum = record[0]['GBSeq_other-seqids'][i].split('|')[1]
                    break
        except:
            print "         {0}: Fail to get gi number".format(accessionNum)
            continue

        try:
            fastaSeq = record[0]['GBSeq_sequence']
        except:
            print "         {0}: Failed to record sequence".format(accessionNum)
            continue

        counter = counter + 1

        #Record phylum and phylum name in lists
        phylumQueue.put(phylum)
        sequencesQueue.put(fastaSeq)
        dataQueue.put([eValue,seqLen,giNum,phylum,accessionNum])

def printWrite(str, file):
    print str,
    file.write(str)
    return
        
def fetch_results(proteinCode,outputseq,hmmBool,fileFlag):
    global br
    global lo_cutoff
    
    print "    fetching results..."
    dataList = [] #stores the tuple [e-Value, sequence length, GI Number, phylum] respectively
    phylumNameList = [] #list containing all names of phylum, no duplicate
    phylumList = [] #list containing all names of phylum, with duplicates
    sequences = [] #stores sequences of current family
    avgLen = 0 #Average length of sequence
    counter = 0
    invalidProteins = 0

    # Have to use beautiful soup to decode since the result page is not
    # that nicely formatted
    try:
        soup = BeautifulSoup(br.response(), "html.parser")
        soup = soup.prettify('utf-8')
    except:
        print "Reading from file Directly"
        soup = br.response().read()

    keyLines = soup.split('Sequences with E-value WORSE than threshold')[0]

    #The link containing the accession number contains a tag title "Show
    #  report for"
    keyLines = keyLines.split('Go to alignment for ')
    del keyLines[0]

    with open(proteinCode + '.txt', 'w') as file:
        printWrite("    {0} BLAST hits found in this family. {1} minutes expected to finish\n".format(len(keyLines), round(0.9 * len(keyLines) / 60), -3),file)
        m = mp.Manager()

        dataQueue = m.Queue()
        phylumQueue = m.Queue()
        sequencesQueue = m.Queue()
        invalidQueue = m.Queue()
        
        jobs = []
        for i in range(0,4):
            p = mp.Process(target=process_data,args=(i,keyLines,dataQueue,phylumQueue,sequencesQueue,invalidQueue,lo_cutoff))
            jobs.append(p)
            p.start()

        for item in jobs:
            item.join()

        maxGiLen = 0
        maxAccessionLen = 0
        while (not(dataQueue.empty()) or not(phylumQueue.empty()) or not(sequencesQueue.empty()) or not (invalidQueue.empty())):
            if(not(dataQueue.empty())):
                data = dataQueue.get()
                if(maxGiLen < len(data[2])): maxGiLen = len(data[2])
                if(maxAccessionLen < len(data[4])): maxAccessionLen = len(data[4])
                dataList.append(data)
            if(not(phylumQueue.empty())):
                phylum = phylumQueue.get()
                phylumList.append(phylum)
                if phylum not in phylumNameList:
                    phylumNameList.append(phylum)
            if(not(sequencesQueue.empty())):
                sequences.append(sequencesQueue.get())
            if(not(invalidQueue.empty())):
                invalidProteins = invalidProteins + invalidQueue.get()

        #Final outputs
        total = 0
        for num in dataList:
            #compute total sequence length
            total = total + num[1]

        #divide for average
        if(len(dataList) > 0):
            avgLen =  total/len(dataList)
        else:
            avgLen = 0

        if(outputseq):
            
            if fileFlag:
                if bool(re.search('/', proteinCode)):
                    tmp = re.search('/([A-Za-z0-9\-\.\_]+)$', proteinCode)
                    proteinCode = tmp.group(1)
                    proteinCode = os.path.splitext(proteinCode)[0]
		    if not os.path.exists('./clustalout'):
                       os.makedirs('clustalout')
            f = open('./clustalout/'+proteinCode+'.faa','w')
            count = 0
            for item in sequences:
                f.write('>{0}\n'.format(dataList[count][4]))
                f.write(item+('\n'))
                count = count+1

        if(outputseq):
            if not os.path.exists('./tmsout'):
                os.makedirs('tmsout')
            f = open('./tmsout/'+proteinCode+'.faa','w')
            count = 0
            for item in sequences:
                f.write('>{0}\n'.format(dataList[count][4]))
                f.write(item+('\n'))
                count = count+1
            f.close()

        tmsOutput = []
        if(hmmBool and outputseq):
            os.system('hmmtop -if=./tmsout/'+proteinCode+'.faa -of=./tmsout/'+proteinCode+'.txt 2> /dev/null')
            hand = open('./tmsout/'+proteinCode+'.txt')
            tmsOutput = re.findall('\s+(IN|OUT)\s+(\d+)',hand.read())
            hand.close()

        tmsGlobalAvg = 0

        #compute standard deviation
        total = 0
        counter = 0
        for item in dataList:
            total = total + (item[1]-avgLen)**2
            if(hmmBool and outputseq):
               tmsGlobalAvg = tmsGlobalAvg + float(tmsOutput[counter][1])
            counter = counter + 1
        if(len(dataList) > 0):
            stddev = math.sqrt(total/len(dataList))

        else:
            stddev = 0
            
        printWrite("\n    {0} BLAST hits were discarded. {1} BLAST hits with E-Value better than cutoff\n".format(invalidProteins,len(keyLines)-invalidProteins),file)
        
        printWrite("\n    \tGI number\tAccession Number\tE-Value\tLength\tNumber of TMS'S\n",file)

        for phylumName in phylumNameList:
            total = phylumList.count(phylumName)
            printWrite('\n    {0} from phylum {1} - {2}%\n'.format(total, phylumName,(float(total)/float(len(dataList)))*100), file)

            #list used for sorting
            phylaData = []
            counter = 0
            tmsavg = 0
            tmstotal = 0
            while counter in range(len(dataList)):
                if (phylumName == dataList[counter][3]):
                    if(hmmBool and outputseq):
                       tms = float(tmsOutput[counter][1])
                       tmsavg = tmsavg + tms
                       tmstotal = tmstotal + 1
                       phylaData.append(dataList[counter]+[tms])
                    else:
                       phylaData.append(dataList[counter])
                counter = counter + 1
            phylaData.sort()
            for item in phylaData:
                if(hmmBool and outputseq):

                   printWrite("        {0:<{1}}\t{2:<{3}}\t\t{4}\t{5}\t{6}\n".format(item[2],abs(maxGiLen),item[4],abs(maxAccessionLen),item[0],item[1],item[5]),file)
                else:
                   printWrite("        {0:<{1}}\t{2:<{3}}\t\t{4}\t{5}\n".format(item[2],abs(maxGiLen),item[4],abs(maxAccessionLen),item[0],item[1]),file)
            if(hmmBool and outputseq):
                phylavg =  tmsavg/tmstotal
	        tphylavg ="%.1f" % phylavg 
	        print '        Phyla TMS Average:{0}'.format(tphylavg)

        #Average Length
	       	tavgLen = "%.1f" % avgLen
	tstddev = "%.1f" % stddev
        tottmsavg = tmsGlobalAvg/len(dataList)
	ttottmsavg = "%.1f" % tottmsavg
        printWrite("\n    Alignment average length: {0} aa".format(tavgLen),file)
        printWrite("\n    Standard Deviation: {0} aa".format(tstddev),file)
        if(outputseq and hmmBool):
          printWrite("\n    TMS Average: {0} TMSs\n\n".format(ttottmsavg),file)
	 # printWrite("\n    TMS Standard Deviation: {0} TMSs\n\n".format(tmsstddev),file)

        fusionLen4x = 4.0 * stddev + avgLen
        fusionLen3x = 3.0 * stddev + avgLen
        fusionLen2x = 2.0 * stddev + avgLen

        #Records index values on list
        fusion4x = []
        fusion3x = []
        fusion2x = []

        counter = 0
        while counter in range(len(dataList)):
            if dataList[counter][1] >= fusionLen4x:
                fusion4x.append(dataList[counter])
            elif dataList[counter][1] >= fusionLen3x:
                fusion3x.append(dataList[counter])
            elif dataList[counter][1] >= fusionLen2x:
                fusion2x.append(dataList[counter])
            counter = counter + 1

        #Sort all fusion proteins by e-value
        fusion4x.sort()
        fusion3x.sort()
        fusion2x.sort()

        if (not fusion4x and not fusion3x and not fusion2x):
            printWrite("    No potential fusion proteins found.\n", file)
        else:
            printWrite("    Potential fusion proteins...\n",file)
            if (fusion4x):
                printWrite("    Listing Accession Numbers of proteins 4 standard deviations greater than the mean:\n",file)
                for item in fusion4x:
                    printWrite("        {0};{1}\n".format(item[4],item[0]),file)
            else:
                printWrite("    No potential fusion proteins found 4 standard deviations from the mean.\n", file)

            if (fusion3x):
                printWrite("    Listing Accession Numbers of proteins 3 standard deviations greater than the mean:\n",file)
                for item in fusion3x:
                    printWrite("        {0};{1}\n".format(item[4], item[0]),file)
            else:
                printWrite("    No potential fusion proteins found 3 standard deviations from the mean.\n", file)
                
            if (fusion2x):
                printWrite("    Listing Accession Numbers of proteins 2 standard deviations greater than the mean:\n",file)
                for item in fusion2x:
                    printWrite("        {0};{1}\n".format(item[4], item[0]),file)
            else:
                printWrite("    No potential fusion proteins found 2 standard deviations from the mean.\n", file)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='This program takes one or more TCDB family codes (i.e. 1.A.1) as input and runs PSI-BLAST for all proteins in each family. The taxanomic classification for each BLAST hit is extractcted and the protein composition of each detected phyla within the query family(-ies) is reported. PhyST can generate the FASTA sequences of top scoring hits, which can be used to run other programs (e.g. multiple alignment, hydropathy analyses, etc). PhyST will also run ClustalW if indicated by the user. For examples go to http://biotools.tcdb.org/barphyst.html'.format(sys.argv[0]))
    parser.add_argument('-e',type=float,metavar='float',help='BLAST E-value. For example, -e 0.0001 or -e  1e-4.',nargs='?',default=0.0001)
    parser.add_argument('-i',type=int,metavar='int',help='Number of PSI-BLAST iterations.',nargs='?',default=2)
    parser.add_argument('-a',type=int,metavar='int',help='Minimum alignment coverage for the query. By default,By default if 50 percent of the protein is not aligned then the match is discarded.',nargs='?',default=50)
    parser.add_argument('-w',type=int,metavar='int',help='The length of the wordsize that initiates a BLAST alignment. Must be 2,3, or 6.',nargs='?',default=3)
    parser.add_argument('-u',help = 'If this option is given the program will run ClustalW for you.',action = 'store_true',default= False)
    parser.add_argument('-f',help = 'Reads files with fasta sequences instead of TCDB ID Numbers.',action = 'store_true',default = False)
    parser.add_argument('family',metavar='Family/Protein',help='A list of space-separated TCDB ID numbers (e.g  1.E.5 1.E.1.1 1.E.1.2). At least 1 TCDB ID must be given.',nargs='+')
   
    args = parser.parse_args()

    if(args.u and (not datatrue)):
        print("Option -u requires option -c. Please add -c to your command  and try again")
        exit() 

    #Set globals
    lo_cutoff = args.a
    eValue = str(args.e)
    if (args.w == 3 or args.w == 2 or args.w == 6):
        word_size = str(args.w)
    else:
        print "Word size is incorrect. Must be either 2,3, or 6.\n"
        exit()

    #begin parsing
    print "Will search {0} families. Query cover cut off: {1}%.".format(len(args.family), args.a)
    familyCount = 0
    for fam in args.family:
        familyCount = familyCount + 1
        retries = 0
        
        while (retries < 2):
            start_time = time.time()
            if(args.f):
              print "\nI am working on file of protein family not on TCDB:".format(familyCount,fam)
            else:
              print "\nI am working on family #{0}: {1}".format(familyCount, fam)
            try:
                blast(fam,args.f)
            except:
                print "    The FASTA  sequence of the family cannot be retrieved in order to conduct PSI-BLAST, please check your netwrok connectivity and tcdb.org to see if family {0} listed.".format(familyCount)
            if(timed_out):
                if(retries < 2):
                    print "    Retrying in 10 seconds..."
                sleep(10)
                timed_out = False
                retries = retries + 1
            else:
                retries = 2
                #do iterations
                print "Iterating %d times:" % args.i
                for i in range(1,args.i+1):
                    print "    performing iteration %d," %(i)
                    iterate()
                fetch_results(fam,datatrue,tmstrue,args.f)
                if(datatrue and args.u):
                     
                     fileID = ""
                     #Handles the option of running BLAST on protein families not listed in TCDB
                     if (args.f):
                         if bool(re.search('/', fam)):
                             tmp = re.search('/([A-Za-z0-9\-\.\_]+)$', fam)
                             fileID = tmp.group(1)
                             trueFile = ntpath.basename(fam)
                         origFile= os.path.splitext(trueFile)[0]
                         os.system("clustalo -i clustalout/" + origFile + ".faa    -o clustalout/" + origFile + ".aln    > /dev/null")
                         print"\n\n\nRan ClustalOmega for family %s and .aln file has been created \n" % fam
         
                     if not fileID:
                         fileID = fam
                             
                     if(args.f == False):
                       os.system("clustalo -i ./clustalout/" + fam + ".faa -o clustalout/" + fam + ".aln  > /dev/null")
                       print"\n\n\nRan ClustalOmega for family %s and .aln file has been created \n" % fam

    close()
            


print "\nProgram finished."
