#!/usr/bin/env python

# ramdisk.py by Vamsee Reddy - This tool will allow you to create a temporary directory that will serve as a ramdisk. 
# Just launch disk=ramdisk.Create(DISKSIZE, PATH_TO_CREATE)
# to unmount run disk.unmount

import os
import subprocess
import shutil

class Create:

    def __init__(self,size,path): 
        self.size=int(size)
        self.ramdir = path
        if self.disk_exists() is False:
            self.create_disk()
        return

    def disk_exists(self):
        if os.path.isdir(self.ramdir) is True:
            return True
        return False

    def create_disk(self):
        cmd = "hdid -nomount ram://%i"%(self.size)
        command=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        id = command.communicate()
        self.diskid = id[0].strip()
        create = 'newfs_hfs %s' %(self.diskid)
        subprocess.Popen(create,shell=True,stdout=subprocess.PIPE).wait()
        os.mkdir(self.ramdir)
        mount = 'mount -t hfs %s %s' %(self.diskid,self.ramdir)
        subprocess.Popen(mount,shell=True,stdout=subprocess.PIPE).wait()
        return

    def unmount(self):
        cmd = "umount %s" %(self.diskid)
        subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).wait()
        shutil.rmtree(self.ramdir)
