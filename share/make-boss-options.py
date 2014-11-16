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
  sharedir=os.path.abspath(os.curdir+"/../share/")
  f.write("#!/bin/tcsh\n")
  f.write("cd "+rundir+"\n")
  f.write("source "+sharedir+"/setup.csh\n")
  f.write("boss.exe "+str(run)+".cfg\n")



def proceed(run, directory, files):
    TEMPLATE_RUN_NUMBER=str(run)
    TARGET_FILE = TEMPLATE_RUN_NUMBER+".cfg"
    r = "run_0+"+str(run)+"_.+\.dst"
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
    configure('../share/template.cfg',TARGET_FILE,TEMPLATE_RUN_NUMBER, TEMPLATE_DST_FILES)
#create qsub files
    configure_pbs_jobs(run)

psi2s2011 = range(25244,25338);
jpsi2011 = range(24937,24979);
tau2011 = range(24984,25244)

for run in tau2011:
  print "Proceeding run ", run
  os.path.walk("data", proceed, run)
