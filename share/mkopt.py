#!/usr/bin/python
import os
import sys
import OptionMaker


from OptionMaker import OptionMaker
from optparse import OptionParser

opt = OptionParser()


opt.add_option("-t", "--template", dest="template_file",
                  help="Template configuration file", metavar="FILE")

opt.add_option("-p", "--prefix", dest="job_prefix",
                  help="Prefix for the pbs job", default="test")

opt.add_option("-j", "--jobs", dest="jobs_number",
                  help="The jobs number")

opt.add_option("-N", "--events", dest="event_number","default=100"
                  help="The jobs number")


(options, args) = opt.parse_args()

JPSIKKROOT_DIR = os.path.abspath(os.environ['JPSIKKROOT'])
TEMPLATE_DIR = os.path.join(JPSIKKROOT_DIR, 'share/template')
TEMPLATE_DECAY_FILE=os.path.join(TEMPLATE_DIR, 'decay.dec')
TEMPLATE_FILE = "selection.cfg"


SelectionMode = False
SimulationMode = False
ReconstructionMode = False

if len(args) < 1:
	print 'Specify the action: "sel",  "sim",  "rec"'
	sys.exit(1)

if args[0] == "selection" or args[0] == "sel":
	SelectionMode = True
	TEMPLATE_FILE = "selection.cfg"
	print "Making selection config files."

if args[0] == "simulation" or args[0] == "sim":
	SimulationMode = True
	TEMPLATE_FILE = "simulation.cfg"
	print "Making simulation config files."

if args[0] == "reconstruction" or args[0] == "rec":
	ReconstructionMode = True
	TEMPLATE_FILE = "reconstruction.cfg"
	print "Making reconstruction config files."


#configure template config file
if not os.path.exists(TEMPLATE_FILE):
	TEMPLATE_FILE = os.path.join(TEMPLATE_DIR, TEMPLATE_FILE)
else:
	TEMPLATE_FILE = os.path.realpath(os.path.abspath(TEMPLATE_FILE))
if not os.path.exists(TEMPLATE_FILE):
	print 'Template file ',  TEMPLATE_FILE,  ' doest exists'
	sys.exit(1)

print "Setup template file: ",  TEMPLATE_FILE

if SimulationMode:
	if len(args) < 3:
		print 'Specify decay file of the simulation and event number'
		sys.exit(1)
	TEMPLATE_DECAY_FILE = args[1]
        TEMPLATE_EVENT_NUMBER=int(args[2])
	if not os.path.exists(TEMPLATE_DECAY_FILE):
		TEMPLATE_DECAY_FILE = os.path.join(TEMPLATE_DIR, TEMPLATE_DECAY_FILE)
		if not os.path.exists(TEMPLATE_DECAY_FILE):
			print 'Decay file ',  TEMPLATE_DECAY_FILE,  ' doest exists'
			sys.exit(1)
	else:
		TEMPLATE_DECAY_FILE = os.path.realpath(os.path.abspath(TEMPLATE_DECAY_FILE))
	print "Use decay file: ",  TEMPLATE_DECAY_FILE
        optMaker = OptionMaker(TEMPLATE_DECAY_FILE, TEMPLATE_EVENT_NUMBER, options.jobs_number)
	TARGET_DIR="."
	if len(args) == 3:
		TARGET_DIR = os.path.realpath(os.path.abspath(args[2]))
        optMaker.make_sim(TEMPLATE_FILE, TARGET_DIR, options.job_prefix)



if SelectionMode:
	if len(args) < 2:
		print 'specify source dir with dst files or dst file'
		sys.exit(1)
	DATA_SOURCE_DIR = args[1]
	if not os.path.exists(DATA_SOURCE_DIR):
			print 'Source data dir or dst file ',  DATA_SOURCE_DIR,  ' does not exists'
	if os.path.isdir(DATA_SOURCE_DIR):
		DATA_SOURCE_DIR = os.path.realpath(DATA_SOURCE_DIR)
	TARGET_DIR="."
	if len(args) == 3:
		TARGET_DIR = os.path.realpath(os.path.abspath(args[2]))
	print "TARGET_DIR =",  TARGET_DIR
	optMaker = OptionMaker(DATA_SOURCE_DIR)
	print "Found ",  len(optMaker.fileList),  " *.dst files"
	optMaker.group()
	print "Found ",  len(optMaker.runMap), " runs,  creating joboptions..."
	optMaker.make(TEMPLATE_FILE,TARGET_DIR,options.job_prefix)

	
#if( args[1] == "simulation" || args[1] == "sim"):
#	print "Making simulation config files"
#
#if( args[1] == "reconstruction" || args[1] == "rec"):
#	print "Making reconstruction config files"

#TEMPLATE_FILE = os.path.join(TEMPLATE_DIR, 'selection.cfg')
#
#print options.selection
#print args
#DATA_SOURCE_DIR = os.path.realpath(os.path.abspath(args[0]))
#
##optMaker.group(".*run_(\d\d\d\d\d\d\d).*.dst")
