#!/usr/bin/python 
#This script create simulation and reconstruction options
import os
import sys
import re

if len(sys.argv)<4:
    print "Usage: make-sim-options.py <decay_file> <output_prefix> <event_number>"
    exit(1)

HOME_DIR = os.environ['HOME']
JPSIKKROOT_DIR = os.environ['JPSIKKROOT']
SHARE_DIR = os.path.join(JPSIKKROOT_DIR, "share")
TEMPLATE_DIR = os.path.join(JPSIKKROOT_DIR, "share/template")

TEMPLATE_SIM_FILE = os.path.joint(TEMPLATE_DIR, "simulation.cfg")
print HOMEDIR, JPSIKKROOT_DIR, TE

DECAY_FILE = os.path.abspath(os.path.join(SHARE_DIR,sys.argv[1]))
PREFIX = sys.argv[2]

RTRAW_FILE = os.path.abspath(PREFIX+".rtraw")
DST_FILE = os.path.abspath(PREFIX+".dst")
ROOT_FILE = os.path.abspath(PREFIX+".root")

