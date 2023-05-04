#!/usr/bin/env python
'''
  phyST
  This program takes a protein family code (i.e. 1.A.1) as input, then
  searches for the alignment sequence of the first protein in that
  family (i.e. 1.A.1.1.1), then psi-blasts this sequence, and finally
  outputs the number of protein homologues from every phylum found as
  well as the average sequence length of these homologues
  NOTE: TO CREATE CLUSTAL OUTPUT (.FLAT FILE) MAKE A FOLDER CALLED
    "clustalout" IN THE SAME DIRECTORY AS THE SCRIPT.
'''
# Written by Hari Krishnan, Larry Chau
# lchau@ucsd.edu - larrymchau@gmail.com
# hkkrishn563@gmail.com - hkkrishn@ucsd.edu

import mechanize
import sys,re,os
import urllib2
import tempfile
from Bio import Entrez
from bs4 import BeautifulSoup
from time import sleep
import cookielib
import math
import multiprocessing as mp
import optparse
import ctypes

#Specify the user's email when using Entrez from Bio package
Entrez.email = "lchau@ucsd.edu"
    
#Globals, working data
interval = 0
br = mechanize.Browser()
tfiles = []
lo_cutoff = 50
eValue = '0.0001'
alignment = ''

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
        seconds = int(time.groups()[0])
        secs = seconds if seconds < 15 else 10
        if(interval >= 15):
            interval = 0
            print '    waited 15 seconds...'
        interval = interval + secs
        sleep(secs)
        return True
    return False

#Code copied from protocol1.py
def loop():
    global br
    if wait():
        br.select_form(nr=0)
        br.submit()
        loop()
    return

def find_optimal():
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

def blast(proteinCode):
    global eValue
    global alignment
    global br
    
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
        link = br.click_link(text_regex = '\d+.\D+.\d.\d.\d', nr= cnt)
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

    #click into the FASTA fornatted sequence, then split texts to
    #  extract the string containing only alignment sequence
    sourcePage = br.open(link)
    keyLines = sourcePage.read().split('<PRE>')[1]
    keyLines = keyLines.split('</PRE>')[0]
    keyLines = keyLines.split('\n')
    del keyLines[0]
    for row in keyLines:
        alignment  = alignment + row
    
    optimEVal = find_optimal()

    print '    Estimate Optimal E-Value (TCDB):',optimEVal,'(e.g. ',float(optimEVal),')'
    print '    Using:',eValue
    eValue = eValue.replace(' ', '')
    
    #Go to NCBI blast page, enter the alignment found above, select
    #  "psi-blast" and "5000" max results", then blast
    br.open('http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome')
    br.select_form(nr=0)
    br.form.set_all_readonly(False)
    br.form['SELECTED_PROG_TYPE'] = 'psiBlast'
    br.form['RUN_PSIBLAST'] = 'on'
    br.form['BLAST_PROGRAMS'] = ['psiBlast']
    br.form['QUERY'] = alignment
    br.form['I_THRESH'] = eValue
    br.form['MAX_NUM_SEQ'] = ['10000']
    br.submit()
    print "    Blasting Off..."
    loop()
    print "    Done.\n"

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
    #results = open('out.html','r').read()
    results = results.replace('<! --','<!--')
    file = open('out.html','w')
    file.write(results)
    
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
                print form.attrs['name']
                print 'id={0}'.format(formcount)
                br.form = form
                break
        formcount=formcount+1
    br.select_form(nr=formcount)
    
    print br.form
    
    '''
        if 'id' in form.attrs:
            if form.attrs['id'] == 'smrtBlastForm':
                print form.attrs['id']
                print 'id={0}'.format(formcount)
                br.form = form
                break
        formcount=formcount+1
        if 'name' in form.attrs:
            if form.attrs['name'] == 'overview0':
                print form.attrs['name']
                print 'id={0}'.format(formcount)
                br.form = form
                break
        formcount=formcount+1
    '''     
    br.form.action='http://blast.ncbi.nlm.nih.gov/Blast.cgi'
    br.submit()

    file = open('out2.html','w')
    file.write(br.response().read())    
 
    loop()
    
    file = open('out3.html','w')
    file.write(br.response().read())
    
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
        
def fetch_results(proteinCode,outputseq):
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
        print "    Trouble soupifying, reading from file directly."
        soup = br.response().read()
        
    keyLines = soup.split('Sequences with E-value WORSE than threshold')[0]

    #The link containing the accession number contains a tag title "Show
    #  report for"
    keyLines = keyLines.split('Go to alignment for ')
    del keyLines[0]

    with open(proteinCode + '.txt', 'w') as file:
        printWrite("    {0} proteins found in this family. {1} minutes expected to finish\n".format(len(keyLines), round(0.9 * len(keyLines) / 60), -3),file)

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
        while (not(dataQueue.empty()) or not(phylumQueue.empty()) or not(sequencesQueue.empty()) or not (invalidQueue.empty())):
            if(not(dataQueue.empty())):
                data = dataQueue.get()
                if(maxGiLen < len(data[2])): maxGiLen = len(data[2])
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
            avgLen = total/len(dataList)
        else:
            avgLen = 0

        if(outputseq):
            if not os.path.exists('./clustalout/'):
                os.mkdir('./clustalout/')
            f = open('./clustalout/'+proteinCode+'.flat','w')
            count = 1
            for item in sequences:
                f.write('>{0}\n'.format(count))
                f.write(item+('\n'))
                count = count+1

        #compute standard deviation
        total = 0
        for item in dataList:
            total = total + (item[1]-avgLen)**2
        if(len(dataList) > 0):
            stddev = math.sqrt(total/len(dataList))
        else:
            stddev = 0
            
        printWrite("    {0} proteins found below the cutoff. Expect {1} proteins.\n".format(invalidProteins,len(keyLines)-invalidProteins),file)
        
        printWrite("\n    \tGI number\tAccession Number\tE-Value\tLength\n",file)
        
        for phylumName in phylumNameList:
            total = phylumList.count(phylumName)
            printWrite('\n    {0} from phylum {1} - {2}%\n'.format(total, phylumName,(float(total)/float(len(dataList)))*100), file)

            #list used for sorting
            phylaData = []
            counter = 0
            while counter in range(len(dataList)):
                if (phylumName == dataList[counter][3]):
                    phylaData.append(dataList[counter])
                counter = counter + 1
            phylaData.sort()
            for item in phylaData:
                printWrite("        {0:<{1}}\t{2}\t\t{3}\t{4}\n".format(item[2],abs(maxGiLen),item[4],item[0],item[1]),file)
        
        #Average Length
        printWrite("\n    Alignment average length: {0} aa".format(avgLen),file)
        printWrite("\n    Standard Deviation: {0} aa \n\n".format(stddev),file)

        fusionLen4x = 4.0 * avgLen
        fusionLen3x = 3.0 * avgLen
        fusionLen2x = 2.0 * avgLen

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

        printWrite("    Potential fusion proteins...\n",file)
        printWrite("    Listing GI Numbers of proteins 4x the avg length:\n",file)
        for item in fusion4x:
            printWrite("        {0};{1}\n".format(item[2],item[0]),file)

        printWrite("    Listing GI Numbers of proteins 3x the avg length:\n",file)
        for item in fusion3x:
            printWrite("        {0};{1}\n".format(item[2], item[0]),file)

        printWrite("    Listing GI Numbers of proteins 2x the avg length:\n",file)
        for item in fusion2x:
            printWrite("        {0};{1}\n".format(item[2], item[0]))

if __name__=='__main__':
    #stores indicies of options
    optionidx = []
    #Default iterate value
    numiterate=0
    #boolean for outputting seq files
    outputseq=0

    #If only 2 arguments detected, output usage
    if len(sys.argv)<2:
        print "Usage:",os.path.basename(__file__),"{[-a[cutoff value %] -i[# iterations] -e[E-Value] -c(write sequences)} [family code...]"
        print "Example:",os.path.basename(__file__),"-a30 -c -e0.001 1.A.1 1.A.2 \n\t{-a30},(-e0.001) and (-c) is optional"
        quit()

    #find command line arg with cutoff
    counter = 0
    for item in sys.argv:
        result = re.match('-a(\d+)',item)
        result2 = re.match('-i(\d+)',item)
        result3 = re.match('-e(\S+)',item)
        result4 = re.match('-c',item)
        if(bool(result)):
            lo_cutoff = int(result.group(0).split('-a')[1])
            optionidx.append(counter)
        if(bool(result2)):
            numiterate = int(result2.group(0).split('-i')[1])
            optionidx.append(counter)
        if(bool(result4)):
            outputseq = 1
            optionidx.append(counter)
        if(bool(result3)):
            eValue = result3.group(0).split('-e')[1]
            optionidx.append(counter)
        counter = counter+1
        
    #begin parsing
    print "Will search {0} families. Query cover cut off: {1}%.".format(len(sys.argv) - len(optionidx) - 1, lo_cutoff)
    cnt = 1
    familyCount = 0
    while cnt in range(len(sys.argv)):
        cnt = cnt + 1
        if (cnt-1 in optionidx):
            continue
        familyCount = familyCount + 1
        proteinCode = sys.argv[cnt-1]
        print "\nI am working on family #{0}: {1}".format(familyCount, proteinCode)
        try:
            blast(proteinCode)
        except:
            print "Error occured on family {0}. So I skip to next family.".format(sys.argv[cnt-1])
            continue
        #do iterations
        print "Iterating %d times:" % numiterate
        for i in range(0,numiterate):
            print "    performing iteration %d," %(i+1)
            iterate()
        fetch_results(proteinCode,outputseq)
        close()

    print "\nProgram finished."
    
