#!/usr/bin/python
import os
import sys
import OptionMaker


from OptionMaker import OptionMaker
from optparse import OptionParser

opt = OptionParser()


opt.add_option("-t", "--template_file", dest="template_file", default="selection.cfg", help="Template configuration file", metavar="FILE")

opt.add_option("-p", "--job_prefix", dest="job_prefix", help="Prefix for the pbs job", default="test")

opt.add_option("-J", "--job_number", dest="job_number",type="int", default=1, help="The jobs number")

opt.add_option("-N", "--event_number", type="int", dest="event_number",default=100, help="Event number per job")

opt.add_option("-r", "--run", dest="runs",default="-8093", help="Number of run in simulation")

opt.add_option("-f", "--force", action="store_true", dest="force", help="Force ")

opt.add_option("-n", "--runs_per_job",type="int", dest="run_number", default=1, help="Number of runs per job")

opt.add_option("-F", "--file_filter",dest="file_filter", default=".*(\d{4,7}).*.dst", help="dst file filter")


(options, args) = opt.parse_args()


optMaker = OptionMaker(options, args)
optMaker.make()

##optMaker.group(".*run_(\d\d\d\d\d\d\d).*.dst")
