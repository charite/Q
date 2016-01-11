#!/usr/bin/env python

import sys
import os.path
import subprocess
import argparse

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

parser = argparse.ArgumentParser(description="Sort and Index BAM files")
parser.add_argument('--output', type=str)
parser.add_argument('input_file')

results, leftovers = parser.parse_known_args()

print results.input_file
#print results.output

if results.output is not None:
 outFilenamePrefix, file_extension = os.path.splitext(results.output)
 outFilenamePrefix = os.path.dirname(results.input_file)+ "/" +outFilenamePrefix;
else:
 outFilenamePrefix, file_extension = os.path.splitext(results.input_file)
 outFilenamePrefix += "_sorted"
#print "output file: " + outFilenamePrefix + ".bam"

inputFile = results.input_file

inFilenamePrefix, inFileExtension = os.path.splitext(results.input_file)
if inFileExtension == ".sam":
 args = ("samtools", "view", "-Sb", results.input_file)
 print "converting to BAM..."
 #print args
 f = open(inFilenamePrefix + ".bam", "w")
 popen = subprocess.Popen(args, stdout=f)
 popen.wait()
 f.close()
 if popen.returncode != 0:
  print "error"
  sys.exit()
 print inFilenamePrefix + ".bam created"
 inputFile = inFilenamePrefix + ".bam"

args = ("samtools", "sort", results.input_file, outFilenamePrefix)
print "sorting..."
#print args
popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output
if popen.returncode != 0:
 print "error"
 sys.exit()
print outFilenamePrefix + ".bam created"
print "indexing..."
args = ("samtools", "index", outFilenamePrefix+".bam")
popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output
if popen.returncode != 0:
 print "error"
 sys.exit()
print outFilenamePrefix+".bai created"