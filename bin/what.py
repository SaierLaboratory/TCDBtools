#!/usr/bin/env python
import hmmtop
import re
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import tempfile
from hashlib import md5
from math import sin,cos,radians,sqrt


class what:

  def __init__(self, seq, window=19, angle=100):
    self.seq = seq
    self.window = int(window)
    self.angle = int(angle)
    self.paths()


  def paths(self):
    index={'G':(-0.400,0.48),'I':(4.500,1.38),'S':(-0.800,-0.18),'Q':(-3.500,-0.85),'E':(-3.500,-0.74),\
          'A':(1.800,0.62),'M':(1.900,0.64),'T':(-0.700,-0.05),'Y':(-1.300,0.26),'H':(-3.200,-0.4),\
          'V':(4.200,1.08),'F':(2.800,1.19),'C':(2.500,0.29),'W':(-0.900,0.81),'K':(-3.900,-1.5),\
          'L':(3.800,1.06),'P':(-1.600,0.12),'N':(-3.500,-0.78),'D':(-3.500,-0.90),'R':(-4.500,-2.53)}
    midpt = (self.window+1)/2
    length = len(self.seq)
    hydro = []
    for i in range(length-self.window+1):
      total = 0
      for j in range(self.window):
        total +=index[self.seq[i+j]][0]
      total = total/self.window
      hydro.append(total)
    self.hydro=hydro

    # Now calculate amphipathicity

    limit = length - self.window + 1
    amphi = []
    for i in range(limit):
      total = 0
      sumsin = 0
      sumcos = 0
      tangle = self.angle

      for j in range(self.window):
        hydro = index[self.seq[i+j]][1]
        sumsin += hydro * sin(radians(tangle))
        sumcos += hydro * cos(radians(tangle))
        tangle = (tangle+self.angle)%360

      sumsin = sumsin * sumsin
      sumcos = sumcos * sumcos
      total = (sqrt(sumsin + sumcos)/self.window)*(5/2) # Correction factor!
      amphi.append(total)
    self.amphi=amphi


  def graph(self, tms = [], title_in = None, save = False):

    #xlim(0, len(self.hydro)) 
    x_data = range(0, len(self.hydro))
    diff=(len(self.seq)-len(self.hydro))/2

    plt.figure()
    plt.axhline(y=0, color='black')
    plt.ylim(-3, 3)
    plt.xlim(right=len(self.seq))
    plt.plot(x_data, self.hydro, linewidth=1.0, label="hydrophobicity", color='b')
    plt.plot(x_data, self.amphi, linewidth=1.0, label="amphipathicity", color='r')

    #plt.legend()
    plt.xlabel("Residue Number")
    plt.ylabel("Value")
    if '/' in title_in:
      plt.title(title_in[title_in.rindex('/')+1:])
    else:
      plt.title(title_in)
    # Highlight TMS regions
    for t in tms:
      plt.axvspan(t[0]-diff,t[1]-diff, facecolor="orange", alpha=0.2)
    width = (0.0265)*len(self.seq) if len(self.seq) > 600 else 15
    if save:
      filename = md5(self.seq).hexdigest()
      fig = matplotlib.pyplot.gcf()
      #fig.frameon = False
      #fig.add_axes([0,0,0,0],frameon=False)
      fig.set_size_inches(width,5.5)
      self.outfile= title_in
      plt.savefig(self.outfile, dpi=80, format="png",bbox_inches='tight', pad_inches=0.003)
    else:
      plt.show()

if __name__=='__main__':

  ht = hmmtop.tools()
  ht.add_library('HT',fastafile.name)
  ht.scan_libraries()
  try:
    tms = ht.results['HT']['what'].values()
  except:
    tms = []
  eh=what(seq,window,angle)
  eh.graph(tms,'Hydropathy & Amphipathicity',True)
  print '<HTML><head><title>Hydropathy & Amphipathicty plot</title></head><body>'
  print "<img src ='%s'><br><br>"%eh.outfile
  print '<font color="Blue">Blue lines denote Hydropathy</font><br>'
  print '<font color="Red">Red lines denote Amphipathicity</font><br>'
  print '<font color="#C5860D">Orange bars mark transmembrane segments as predicted by HMMTOP</font>'
  print '</body></HTML>'
