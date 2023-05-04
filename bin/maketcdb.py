#!/usr/bin/env python
import urllib2
import urllib
import subprocess
TCDB_PATH = 'tcdb'
urllib.urlretrieve ('http://www.tcdb.org/api.php?tcid=all', TCDB_PATH)
subprocess.Popen("makeblastdb -dbtype prot -in " + TCDB_PATH + " -logfile " + TCDB_PATH \
    + ".log", shell=True)
print "done"
