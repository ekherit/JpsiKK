#!/usr/bin/python 
#This script create massive jobOptions files for all dst files
#in data directory grouping them per single run
#
# Usage: make-boss-options.py <data_source_dir> <template_dir> [<target_dir>] [<prefix>]
#
#
import os
import sys
import re

if len(sys.argv)<3:
    print "Usage: make-batch.py <data_source_dir> <template_dir> [<target_dir>] [<prefix>]"
    exit(1)

TEMPLATE_FILE="selection.cfg"
JOB_PREFIX = "JpsiKK"

DATA_SOURCE_DIR = os.path.abspath(sys.argv[1])
TEMPLATE_DIR = os.path.abspath(sys.argv[2])
TARGET_DIR = os.path.abspath(os.curdir)

if len(sys.argv)>=4:
    TARGET_DIR = os.path.abspath(sys.argv[3])

if os.path.exists(DATA_SOURCE_DIR)==False:
    print "Data dir: ", DATA_SOURCE_DIR, " does not exist"
    exit(1)

if os.path.exists(TEMPLATE_DIR)==False:
    print "Template dir: ", TEMPLATE_DIR, " does not exists"
    exit(1)

if os.path.exists(TARGET_DIR)==False:
    os.mkdir(TARGET_DIR)


if len(sys.argv)>=5:
    JOB_PREFIX = sys.argv[4]


#create all file list in directory
def create_file_list(filelist, directory, files):
    for file in files:
        filelist += [os.path.join(directory,file)]
    filelist.sort()

#create dst filtered file list
def create_dst_file_list(filelist, directory, files):
    r = re.compile("run_.+\.dst")
    for file in files:
        if re.match(r,file):
            filelist += [os.path.join(directory,file)]
    filelist.sort()

#filter dst files. Now it is not used
def dst_filter(files):
    r = ".*run_.+\.dst"
    filtered_file_list = []
    for file in files:
        if re.match(r,file):
            filtered_file_list += [file]
    return filtered_file_list


def create_run_dict(files):
    RunMap = {}
    r = re.compile(".*run_(\d\d\d\d\d\d\d).+.dst")
    for file in files:
        m = re.match(r,file)
        if m:
            run = m.group(1)
            if run in RunMap:
                RunMap[run]+=[file]
            else:
                RunMap[run]=[file]
    return RunMap

#create comma separated and qouted list of files
def make_files_string(files):
    files_string=''
    comma=''
    for file in files:
        files_string=files_string+comma+'"'+file+'"'
        comma=',\n  '
    return files_string

#create job option file from template
def create_job_option_file(run, files, templateDir,targetDir):
    #define input and output files in joboptions
    TEMPLATE_INPUT_FILE = make_files_string(files)
    TEMPLATE_OUTPUT_FILE = os.path.join(targetDir,"%s-%s.root" % (JOB_PREFIX,run))
    source_file_name = os.path.join(templateDir,TEMPLATE_FILE)
    #define the name of cfg file
    target_file_name = os.path.join(targetDir, ("%s-%s.cfg" % (JOB_PREFIX, run)))
    print "  ",target_file_name
    source_file = open(source_file_name, 'r')
    target_file = open(target_file_name, 'w')
    for line in source_file:
        line = re.sub("TEMPLATE_INPUT_FILE", TEMPLATE_INPUT_FILE, line)
        line = re.sub("TEMPLATE_OUTPUT_FILE",TEMPLATE_OUTPUT_FILE, line)
        target_file.write(line)
    source_file.close()
    target_file.close()


print "DATA_SOURCE_DIR = ", DATA_SOURCE_DIR
print "TEMPLATE_DIR    = ", TEMPLATE_DIR
print "TARGET_DIR      = ", TARGET_DIR

#create all filelist from source dir
filelist = []
os.path.walk(DATA_SOURCE_DIR, create_dst_file_list, filelist)
#create per run file list
Runs = create_run_dict(filelist)

#for every run create configuration file
print "Creating config files: "
for run, files in Runs.items():
    create_job_option_file(run,files, TEMPLATE_DIR, TARGET_DIR)

