#!/usr/bin/python
import os
import string
import sys
import fileinput
import re


import os

def configure(source_file_name,  target_file_name,  TEMPLATE_RUN_NUMBER, TEMPLATE_DST_FILES):
#Open source template files with simulation configuration
	source_file = open(source_file_name, 'r')
	target_file = open(target_file_name, 'w')
	for line in source_file:
		line = re.sub("TEMPLATE_DST_FILES", TEMPLATE_DST_FILES, line)
		line = re.sub("TEMPLATE_RUN_NUMBER",TEMPLATE_RUN_NUMBER, line)
		target_file.write(line)
	source_file.close()
	target_file.close()

def configure_pbs_jobs(run):
  filename = str(run)+".sh"
  f = open(filename, 'w')
  rundir=os.path.abspath(os.curdir)
  f.write("#!/bin/tcsh\n")
  f.write("cd "+rundir+"\n")
  f.write("source /ihepbatch/bes/nikolaev/bin/boss664\n")
  RunBoss = "boss.exe JpsiKK-%07d.cfg" % run
  f.write(RunBoss)



def proceed(run, directory, files):
    TEMPLATE_RUN_NUMBER = "JpsiKK-%07d" % run
    TARGET_FILE = TEMPLATE_RUN_NUMBER+".cfg" 
    r = "run_%07d_.+\.dst" %run
    flist = []
    TEMPLATE_DST_FILES=''
    for f in files:
      if re.match(r,f):
        name = os.path.join(directory, f);
        flist.append(name)
        if TEMPLATE_DST_FILES=='': comma=''
        else: comma=',\n'
        TEMPLATE_DST_FILES=TEMPLATE_DST_FILES+comma+'"'+name+'"'
    if TEMPLATE_DST_FILES=='': return
    print TEMPLATE_DST_FILES
    configure('/afs/ihep.ac.cn/users/n/nikolaev/batch/6.6.4.p03/JpsiKK/JpsiKK-00-00-01/share/template.cfg',TARGET_FILE,TEMPLATE_RUN_NUMBER, TEMPLATE_DST_FILES)
#create qsub files
    configure_pbs_jobs(run)

psip2009 = range(0,9999999);
test = range(26996,27000);

for run in test :
  print "Proceeding run ", run
  os.path.walk("data", proceed, run)
