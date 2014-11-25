#!/usr/bin/python
import os
import sys
import OptionMaker


from optparse import OptionParser

opt = OptionParser()
opt.add_option("-t", "--template", dest="selection.cfg",
                  help="Template configuration file", metavar="FILE")

opt.add_option("-p", "--prefix", dest="JPKK",
                  help="Prefix for the pbs job", metavar="FILE")

opt.add_option("-s", "--selection",
                  action="store_false", dest="selection", default=True,
                  help="Use selection template")

(options, args) = opt.parse_args()


JPSIKKROOT_DIR = os.path.abspath(os.environ['JPSIKKROOT'])
TEMPLATE_DIR = os.path.join(JPSIKKROOT_DIR, 'share/template')
TEMPLATE_FILE = os.path.join(TEMPLATE_DIR, 'selection.cfg')

print options.selection
print args
DATA_SOURCE_DIR = os.path.abspath(args[0])

optMaker = mkbopt.OptionMaker(DATA_SOURCE_DIR)
#optMaker.group(".*run_(\d\d\d\d\d\d\d).*.dst")
optMaker.group()
optMaker.make(TEMPLATE_FILE,"tmp","test")
