#!/usr/bin/env python3
# Download CDD image from NCBI
import urllib.request, urllib.parse, urllib.error
import re

def fetch(seq,out,width=970):
    packet={'seqinput':seq,'db':'cdd','value':'0.01','filter':'T','maxhits':'500','mode':'rep'}
    try:
        handle = urllib.request.urlopen('http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi',data=urllib.parse.urlencode(packet))
    except:
        return False
    pattern = re.compile(r'var strDataCacheHandle = "(.+?)";')
    session = pattern.search(handle.read())
    session = session.groups()[0]
    url = 'http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?dhandle=%s&show_feat=true&mode=rep&gwidth=%i&output=graph'%(session,width)
    for i in range(3):
        try:
            urllib.request.urlretrieve(url,out)
            return True
        except:
            pass
    return False

if __name__=='__main__':

    seq='MHYQNKMIMDFQLLQEENSKQQTLLKQQISPEIQKNSLQTYNSSEDDENVQDENWDHIYN\
    VPNMGNITKKVLIKLITASIIAFLFLIAEVTGGILAASLAILSDAAHMFSDISGFFISIF\
    SVWIGTKPASTQLSYGYHRSEVIGAMASIFIIWGLTILLLYEATHRIIKQEKVEEPLYML\
    ITAGFGLFCNIIMAKVLHSAPGHSHHGCSHGHSHDHDHDHEHEKDHEHKKNHEHNHNHDH\
    NHNHKHKHKHNHDHSHNSDHNHDHEHNHSHNHNHEHNQGHNNNHGHEHNNEKQKNQIINN\
    KKSGSNCSEFSAKQYGANKQDIAVSLQPIDQNQNQQIQQKNKRKNLSQISAKDNYNLRAA\
    MIHVIGDIIQSIGVLIAALLIYFLDEKTKYIHLADPICTYLFSVLVLFTTFPVAKECIKV\
    LMEGTPTDINIKQFEAELNAIKDIEEIHDLHIWSLSKGKPSLSCHIFCKDNPKEVLKKAT\
    RLCRYYGIYHTTIQVEDYQDKGSKNYIKCDHNIHH'

    fetch(seq,'out.png')
