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
    JPSIKKROOT_DIR       = os.path.abspath(os.environ['JPSIKKROOT'])
    TEMPLATE_DIR         = os.path.join(JPSIKKROOT_DIR, 'share/template')
    dataDir = ""   #folder where data files will be searched
    templateFile=""# the path to the template file
    decayFile=""   #the decay file
    targetDir="."   #the folder where options file will be located
    jobPrefix=""   #prefix of the job
    fileFilter=""  #dst regex filter
    runFilter=""   #regexp run grouping
    eventNumber=-1 #Event number per job
    runNumber=1    #Run number per job
    runs=""
    jobNumber=1    #number of jobs
    Seed=1         #the random seed
    fileList=[]
    runMap={}
    SelectionMode = True
    SimulationMode = False
    Reconstruction = False


    def __init__(self, options, args):
        self.decayFile = os.path.realpath(os.path.abspath(args[1]))
        self.dataDir   = os.path.realpath(os.path.abspath(args[1]))
        self.eventNumber=options.event_number
        self.jobNumber=int(options.job_number)
        self.jobPrefix=options.job_prefix
        self.runs=options.runs

        if len(args) < 1:
            print 'Specify the action: "sel",  "sim",  "rec"'
            sys.exit(1)

        if(args[0] == "selection" or args[0]=="sel"):
            self.SelectionMode=True
            self.fileFilter=file_filter
            self.fileList = filter_file_list(create_file_list(self.dataDir), self.fileFilter)
            self.templateFile = "selection.cfg"

        if args[0] == "simulation" or args[0] == "sim":
            self.SimulationMode = True
            self.templateFile = "simulation.cfg"
            print "Making simulation config files."

        if args[0] == "reconstruction" or args[0] == "rec":
            self.ReconstructionMode = True
            self.templateFile = "reconstruction.cfg"
            print "Making reconstruction config files."


        if not os.path.exists(self.templateFile):
            self.templateFile = os.path.join(self.TEMPLATE_DIR, self.templateFile)
            print self.templateFile
        else:
            self.templateFile = os.path.realpath(os.path.abspath(self.templateFile))
            print self.templateFile
        if not os.path.exists(self.templateFile):
            print 'Template file ',  options.template_file,  ' doest exists'
            sys.exit(1)
        print "Setup template file: ",  self.templateFile
	if len(args) >= 2:
		self.targetDir = os.path.realpath(os.path.abspath(args[2]))
        if not os.path.exists(self.targetDir):
            os.mkdir(self.targetDir)
        else:
            if not options.force:
                print "Target dir ",  self.targetDir,  " already exists, exiting"
                sys.exit(1)



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

    def make(self):
        if self.SimulationMode:
            self.make_sim()




    def make_sel(self, template_file, target_dir, job_prefix):
        #self.jobPrefix = job_prefix
        #self.setup_template_file(template_file, target_dir)

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

    def make_sim(self):
        for job in range(0, self.jobNumber):
            #define input and output files in joboptions
            TemplateOutputFile = os.path.abspath(os.path.join(self.targetDir,"%s-%07d.rtraw" % (self.jobPrefix,self.jobNumber)))
            TemplateInputFile =  self.decayFile 
            #define the name of cfg file
            target_file_name = os.path.join(self.targetDir,"%s-sim-%07d.cfg" % (self.jobPrefix, job))
            target_file_name_rec = os.path.join(self.targetDir,"%s-rec-%07d.cfg" % (self.jobPrefix, job))
            source_file = open(self.templateFile, 'r')
            target_file = open(target_file_name, 'w')
            TemplateRandomSeed = str(random.randint(0,2**32))
            for line in source_file:
                line = re.sub("TEMPLATE_INPUT_FILE", TemplateInputFile, line)
                line = re.sub("TEMPLATE_OUTPUT_FILE",TemplateOutputFile, line)
                line = re.sub("TEMPLATE_RANDOM_SEED",TemplateRandomSeed, line)
                line = re.sub("TEMPLATE_EVENT_NUMBER",str(self.eventNumber), line)
                line = re.sub("TEMPLATE_RUN_NUMBER",self.runs, line)
                target_file.write(line)
            source_file.close()
            target_file.close()
