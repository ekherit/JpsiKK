#Module for configuring boss options
import os
import sys
import re
import random

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

#for every run create file list
def create_run_dict(files,reg):
    RunMap = {}
    r = re.compile(reg)
    for file in files:
        m = re.match(r,file)
        if m:
            #g = int(m.group(1))
            #run = "%07d" % int(m.group(1))
            run = int(m.group(1))
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



class OptionMaker:
    dataDir = ""   #folder where data files will be searched
    templateFile=""# the path to the template file
    targetDir=""   #the folder where options file will be located
    jobPrefix=""   #prefix of the job
    fileFilter=""  #dst regex filter
    runFilter=""   #regexp run grouping
    eventNumber=-1 #Event number per job
    runNumber=1    #Run number per job
    jobNumber=1    #number of jobs
    Seed=1         #the random seed
    fileList=[]
    runMap={}
    decayFile=""

    def __init__(self, decay_file, jobs_number):
        self.jobNumber=jobs_number
        self.decayFile = os.path.realpath(os.path.abspath(decay_file))


    def __init__(self, data_dir, file_filter=".+.dst$"):
        self.dataDir = os.path.abspath(data_dir)
        self.fileFilter=file_filter
        self.fileList = filter_file_list(create_file_list(self.dataDir), self.fileFilter)

    #def group(self, run_filter=".*run_(\d\d\d\d\d\d\d).*.dst", run_number=1):
    def group(self, run_filter=".*run.*(\d{4,7})[^\d]+.*.dst", run_number=1):
        self.runFilter=run_filter
        self.runNumber = run_number
        self.runMap = create_run_dict(self.fileList, run_filter)

    def setup_template_file(self, template_file, target_dir):
        self.templateFile  = template_file
        if not os.path.exists(self.templateFile):
            print "Template file ", self.templateFile, " does not exists"
        self.targetDir  = target_dir
        if not os.path.exists(self.targetDir):
            os.mkdir(self.targetDir)
        else:
            print "Target dir ",  self.targetDir,  " already exists, exiting"
            sys.exit(1)


    def make(self, template_file, target_dir, job_prefix):
        self.jobPrefix = job_prefix
        setup_template_file(template_file, target_dir)

        for run, files in self.runMap.items():
            #define input and output files in joboptions
            TemplateOutputFile = os.path.abspath(os.path.join(self.targetDir,"%s-%07d.root" % (self.jobPrefix,run)))
            TemplateInputFile = make_files_string(files)
            #define the name of cfg file
            target_file_name = os.path.join(self.targetDir,"%s-%07d.cfg" % (self.jobPrefix, run))
            source_file = open(self.templateFile, 'r')
            target_file = open(target_file_name, 'w')
            TemplateRandomSeed = str(random.randint(0,2**32))
            for line in source_file:
                line = re.sub("TEMPLATE_INPUT_FILE", TemplateInputFile, line)
                line = re.sub("TEMPLATE_OUTPUT_FILE",TemplateOutputFile, line)
                line = re.sub("TEMPLATE_RANDOM_SEED",TemplateRandomSeed, line)
                line = re.sub("TEMPLATE_EVENT_NUMBER",self.eventNumber, line)
                line = re.sub("TEMPLATE_RUN_NUMBER",self.runNumber, line)
                target_file.write(line)
            source_file.close()
            target_file.close()

    def make_sim(self, template_file, target_dir, job_prefix):
        self.jobPrefix = job_prefix
        setup_template_file(template_file, target_dir)
        for job in range(0, self.jobNumber):
            #define input and output files in joboptions
            TemplateOutputFile = os.path.abspath(os.path.join(self.targetDir,"%s-%07d.rtraw" % (self.jobPrefix,self.jobNumber)))
            TemplateInputFile =  self.decayFile 
            #define the name of cfg file
            target_file_name = os.path.join(self.targetDir,"%s-sim-%07d.cfg" % (self.jobPrefix, self.jobNumber))
            target_file_name_rec = os.path.join(self.targetDir,"%s-rec-%07d.cfg" % (self.jobPrefix, self.jobNumber))
            source_file = open(self.templateFile, 'r')
            target_file = open(target_file_name, 'w')
            TemplateRandomSeed = str(random.randint(0,2**32))
            for line in source_file:
                line = re.sub("TEMPLATE_INPUT_FILE", TemplateInputFile, line)
                line = re.sub("TEMPLATE_OUTPUT_FILE",TemplateOutputFile, line)
                line = re.sub("TEMPLATE_RANDOM_SEED",TemplateRandomSeed, line)
                line = re.sub("TEMPLATE_EVENT_NUMBER",self.eventNumber, line)
                line = re.sub("TEMPLATE_RUN_NUMBER",-run, line)
                target_file.write(line)
            source_file.close()
            target_file.close()
