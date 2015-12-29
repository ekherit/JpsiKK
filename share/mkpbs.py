#!/usr/bin/python

import os
#import stat
import sys
import re
import random
from optparse import OptionParser
from stat import *

#create all file list (used for path.walk
def proceed_create_file_list(filelist, directory, files):
    for file in files:
        filelist += [os.path.join(directory,file)]
    filelist.sort()

#create all file list in directory
def create_file_list(directory):
    files = []
    os.path.walk(directory, proceed_create_file_list, files)
    return files;

def filter_file_list(files, reg):
    r = re.compile(reg)
    filtered_file_list = []
    for file in files:
        if re.match(r,file):
            filtered_file_list += [file]
    return filtered_file_list


def group_files(file_list, n) :
    big_file_list=[]
    count=0
    tmp_list=[]
    for f in file_list:
        if count <n:
            tmp_list.append(f)
        else:
            big_file_list.append(tmp_list)
            tmp_list=[]
            count = 0
            tmp_list.append(f)
        count=count+1
    big_file_list.append(tmp_list)
    return big_file_list


options = OptionParser()
options.add_option("-p", "--job_prefix", dest="job_prefix", help="Prefix for the pbs job", default="test")
options.add_option("-n", "--runs_per_job",type="int", dest="run_number", default=1, help="Number of runs per job")
options.add_option("-q", "--queue",dest="queue", default="besq", help="Queue name")
(opt, args) = options.parse_args()

dir=args[0]
absdir = os.path.abspath(dir)
files = filter_file_list(create_file_list(dir),".+.cfg$")
files.sort()
groups = group_files(files,opt.run_number)

pbs_file_list = []
i=0
for flist in groups:
    tcsh_file_name = "%s-%04d.tcsh" % (opt.job_prefix, i) 
    pbs_name = "%s/%s" % (dir , tcsh_file_name )
    pbs_file = open(pbs_name, 'w')
    s="""#!/bin/tcsh
#PBS -N """ + tcsh_file_name + """
#PBS -o """ + tcsh_file_name + """.log
#PBS -j oe
##PBS -q besq
source /ihepbatch/bes/nikolaev/bin/boss664
cd """ + absdir + """
"""
    for f in flist:
        log = os.path.abspath(os.path.splitext(f)[0]+".log")
        s = s+ "boss.exe "+os.path.abspath(f)+" >& " + log + "\n"
    pbs_file.write(s)
    pbs_file_list.append(pbs_name)
    os.chmod(pbs_name, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)
    i=i+1

submit_file_name = dir+"/"+"submit.sh"
submit_file = open(submit_file_name, 'w')
s = "#!/bin/bash\n"
for pbs_file  in pbs_file_list:
    s = s + "qsub -q " + opt.queue + " " +  os.path.abspath(pbs_file) + "\n"
submit_file.write(s)
os.chmod(submit_file_name, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)

