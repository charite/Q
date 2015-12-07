#!/usr/bin/env python
# Author: Benjamin Menkuec
# Copyright 2015 Benjamin Menkuec
# License: LGPL

import sys
import os.path
import subprocess
import argparse
import platform
import multiprocessing

# this function calls a subroutine that will call samtools to 
# sort the bam file and then build an index of it
def indexBamFile(filename, script_path):
	args = ("python", script_path + "/bam_indexer.py", filename)
	print args
	print "Creating indexed bam file..."
	popen = subprocess.Popen(args, stdout=subprocess.PIPE)
	popen.wait()
	output = popen.stdout.read()
	if results.verbose == True:
		print output
	if popen.returncode != 0:
		print "error"
		sys.exit()
	return;

script_path = os.path.dirname(os.path.realpath(__file__))
dataDir = os.getcwd() + "/"

parser = argparse.ArgumentParser(description="Preprocess fastq files and do mapping")
parser.add_argument('--adapters', type=str, default = "/../data/adapters.fa")
parser.add_argument('--barcodes', type=str, default = "/../data/barcodes.fa")
parser.add_argument('--flexcat_er', type=str, default = "0.2")
parser.add_argument('--flexcat_ol', type=str, default = "4")
parser.add_argument('--flexcat_fm', type=str, default = "19")
parser.add_argument('--flexcat_ml', type=str, default = "19")
parser.add_argument('--flexcat_oh', type=str, default = "0")
parser.add_argument('--flexcat_times', type=str, default = "1")

parser.add_argument('--exo', action='store_true')
parser.add_argument('--clean', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--bowtie_location', nargs='?', default = "")
parser.add_argument('--num_threads', nargs='?', default = 
	str(multiprocessing.cpu_count()))
parser.add_argument('input_file')
parser.add_argument('--output_dir', type=str)
parser.add_argument('--filter_chromosomes', type=str, default="(.*)[H|U|M|_]+(.*)")
parser.add_argument('--random_split', action='store_true')
parser.add_argument('genome', type=str)

results, leftovers = parser.parse_known_args()

flexcatAdapterFilename = (os.path.dirname(os.path.realpath(__file__)) + 
	results.adapters)
flexcatBarcodeFilename = (os.path.dirname(os.path.realpath(__file__)) + 
	results.barcodes)

print "Reads: " + results.input_file
print "Genome: " + results.genome

genomeFilename = results.genome
bowtieLocation = results.bowtie_location
if len(bowtieLocation) > 0:
    bowtieLocation = bowtieLocation + "/"

if(platform.system() == "Windows" and results.bowtie_location == ""):
 print "Bowtie location is required under windows"
 sys.exit()

inputFile = os.path.abspath(results.input_file)

temp, inFileExtension = inputFile.split(os.extsep, 1)
# The output file of flexcat can not be a zipped fastq,
# because bowtie can not handle it.
# Therefore, remove .gz ending if its a fastq.gz file
inFileExtension, temp = os.path.splitext(inFileExtension) 
inFileExtension = "." + inFileExtension
inFilenamePrefixWithoutPath = os.path.basename(inputFile)
inFilenamePrefixWithoutPath, temp = inFilenamePrefixWithoutPath.split(os.extsep, 1)

# set the output filename of flexcat
if results.output_dir is not None:
    outputDir = os.path.abspath(results.output_dir) 
    head, tail = os.path.split(results.output_dir)
    if len(tail) > 0: 
        inFilenamePrefixWithoutPath = tail;
    outputDir = outputDir + "/"
else:
    outputDir = inFilenamePrefixWithoutPath
flexcatOutputFilename = (outputDir + "/" + inFilenamePrefixWithoutPath + 
	inFileExtension)

# if the exo option is set, there wont be any fixed barcode matching.
# therefore the output file of flexcat does not have a postfix
if results.exo:
 bowtieInputFilename = (outputDir + "/" + inFilenamePrefixWithoutPath + 
	inFileExtension)
else:
 bowtieInputFilename = (outputDir + "/" + inFilenamePrefixWithoutPath + 
	"_matched_barcode" + inFileExtension)
bowtieOutputFilename = outputDir + "/" + inFilenamePrefixWithoutPath + ".sam"

# platform independant way of calling flexcat
flexcat_path = script_path+"/../bin/flexcat"
if(platform.system() == "Windows"):
	flexcat_path += ".exe"
	
if results.exo:
 args = (flexcat_path, results.input_file, "-tt", "-t", "-ss", "-st", "-app", 
	"-tnum", results.num_threads, "-times", results.flexcat_times, "-er", 
	results.flexcat_er, "-ol", results.flexcat_ol, "-oh", results.flexcat_oh, 
	"-fm", results.flexcat_fm, "-ml", results.flexcat_ml, "-a", 
	flexcatAdapterFilename,"-o", flexcatOutputFilename)
else:
 args = (flexcat_path, results.input_file, "-tl", "5", "-tt", "-t", "-ss", "-st", 
	"-app", "-tnum", results.num_threads, "-times", results.flexcat_times, "-er", 
	results.flexcat_er, "-ol", results.flexcat_ol, "-oh", results.flexcat_oh, 
	"-fm", results.flexcat_fm, "-ml", results.flexcat_ml, "-b", 
	flexcatBarcodeFilename, "-a", flexcatAdapterFilename,"-o", flexcatOutputFilename)
if not os.path.exists(outputDir):
 os.makedirs(outputDir)
if (os.path.isfile(bowtieInputFilename) == False or results.overwrite == True):
    print "Filtering pre-mapping barcodes and trimming adapters..."
    if results.verbose == True:
        popen = subprocess.Popen(args + tuple(leftovers))
    else:
        popen = subprocess.Popen(args + tuple(leftovers), stdout=subprocess.PIPE)
    popen.wait()
    if popen.returncode != 0:
     print "error"
     sys.exit()
    if results.verbose == True:
        print flexcatOutputFilename + " created"

head, tail = os.path.split(genomeFilename)
genomeIndex, file_extension = os.path.splitext(tail)
 
# check if bowtie index already exists
# if yes, skip building the index
genomeIndexFile = os.path.dirname(genomeFilename) + "/" + genomeIndex;
if (os.path.isfile(genomeIndexFile + ".1.ebwt") == False):
 if(platform.system() == "Linux" or platform.system() == "Linux2"):
  args = (bowtieLocation + "bowtie-build", "-o", "1", genomeFilename, 
	genomeIndexFile)
 else:
  args = ("python",  bowtieLocation + "bowtie-build", "-o", "1", 
	genomeFilename, genomeIndexFile)
 popen = subprocess.Popen(args, stdout=subprocess.PIPE)
 popen.wait()
 output = popen.stdout.read()
 if results.verbose == True:
    print output
 if popen.returncode != 0:
  print "error"
  sys.exit() 

# call bowtie if mapped file does not exist or overwrite is true
if (os.path.isfile(bowtieOutputFilename) == False or results.overwrite == True):
    args = (bowtieLocation + "bowtie", "-S", "-p", results.num_threads, 
		"--chunkmbs", "512", "-k", "1", "-m", "1", "-v", "2", "--strata", "--best", 
		genomeIndexFile, bowtieInputFilename, bowtieOutputFilename)
    if(platform.system() != "Linux" and platform.system() != "Linux2"):
     args = ("python",) + args
    print "Mapping reads..."
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    if results.verbose == True:
        print output
    if popen.returncode != 0:
     print "error"
     sys.exit()
 
# nexus-pre
nexusOutputFilename = (outputDir + "/" + inFilenamePrefixWithoutPath + 
	"_filtered.bam")
if results.random_split == True:
	nexusOutputFilenameSplit1 = (outputDir + "/" + inFilenamePrefixWithoutPath + 
	"_filtered_split1.bam")
	nexusOutputFilenameSplit2 = (outputDir + "/" + inFilenamePrefixWithoutPath + 
	"_filtered_split2.bam")

# platform independant way of calling nexus-pre
nexcat_path = script_path+"/../bin/nexcat"
if(platform.system() == "Windows"):
	nexcat_path += ".exe"
	
# call nexcat if overwrite is true or output files dont exist	
if (os.path.isfile(nexusOutputFilename) == False or 
	(results.random_split == True and 
		(os.path.isfile(nexusOutputFilenameSplit1) == False or 
			os.path.isfile(nexusOutputFilenameSplit2) == False)) or
		results.overwrite == True):
    args = (nexcat_path, bowtieOutputFilename,  "-fc", results.filter_chromosomes)
    if results.random_split == True:
		args += ("-rs",)
    print "Filtering post-mapping barcodes..."
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    if results.verbose == True:
        print output
    if popen.returncode != 0:
        print "error"
        sys.exit()

# special option for binding characteristic analysis
indexBamFile(nexusOutputFilename, script_path)			
if results.random_split == True:
	indexBamFile(nexusOutputFilenameSplit2, script_path)		
	indexBamFile(nexusOutputFilenameSplit1, script_path)		
	
# cleanup
if results.clean:
    print "deleting intermediate files..."
    os.remove(bowtieOutputFilename)
    os.remove(nexusOutputFilename)
    if results.exo:
        os.remove(flexcatOutputFilename)
    else:
        os.remove(outputDir + "/" + inFilenamePrefixWithoutPath + 
			"_matched_barcode" + inFileExtension)
        os.remove(outputDir + "/" + inFilenamePrefixWithoutPath + 
			"_unidentified" + inFileExtension)
