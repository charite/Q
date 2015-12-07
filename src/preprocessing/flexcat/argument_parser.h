// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================
#pragma once

#include <string>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#include "flexi_program.h"

#include "adapter_trimming.h"
#include "demultiplex.h"
#include "general_processing.h"

// Classes

class ArgumentParserBuilder
{
public:
    virtual seqan::ArgumentParser build() = 0;

protected:
    void addGeneralOptions(seqan::ArgumentParser & parser, const FlexiProgram flexiProgram);
    void addFilteringOptions(seqan::ArgumentParser & parser);
    void addDemultiplexingOptions(seqan::ArgumentParser & parser);
    void addAdapterTrimmingOptions(seqan::ArgumentParser & parser, const FlexiProgram flexiProgram);
    void addReadTrimmingOptions(seqan::ArgumentParser & parser, const FlexiProgram flexiProgram);
};

class FilteringParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build() override
    {
        seqan::ArgumentParser parser("sflexFilter");
        addHeader(parser);
        addGeneralOptions(parser, FlexiProgram::FILTERING);
        addFilteringOptions(parser);
        return parser;
    }
private:
    void addHeader(seqan::ArgumentParser & parser);
};

class AdapterRemovalParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build() override
    {
        seqan::ArgumentParser parser("sflexAR");
        addHeader(parser);
        addGeneralOptions(parser, FlexiProgram::ADAPTER_REMOVAL);
        addAdapterTrimmingOptions(parser, FlexiProgram::ADAPTER_REMOVAL);
        return parser;

    }
private:
    void addHeader(seqan::ArgumentParser & parser);
};

class DemultiplexingParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build() override
    {
        seqan::ArgumentParser parser("sflexDMulti");
        addHeader(parser);
        addGeneralOptions(parser, FlexiProgram::DEMULTIPLEXING);
        addDemultiplexingOptions(parser);
        return parser;
    }

private:
    void addHeader(seqan::ArgumentParser & parser);
};

class QualityControlParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build() override
    {
        seqan::ArgumentParser parser("sflexQC");
        addHeader(parser);
        addGeneralOptions(parser, FlexiProgram::QUALITY_CONTROL);
        addReadTrimmingOptions(parser, FlexiProgram::QUALITY_CONTROL);
        return parser;
    }

private:
    void addHeader(seqan::ArgumentParser & parser);
};

class AllStepsParserBuilder : public ArgumentParserBuilder
{
public:
    seqan::ArgumentParser build() override
    {
        seqan::ArgumentParser parser("Flexcat");
        addHeader(parser);
        addGeneralOptions(parser, FlexiProgram::ALL_STEPS);
        addFilteringOptions(parser);
        addDemultiplexingOptions(parser);
        addAdapterTrimmingOptions(parser, FlexiProgram::ALL_STEPS);
        addReadTrimmingOptions(parser, FlexiProgram::ALL_STEPS);
        return parser;
    }

private:
    void addHeader(seqan::ArgumentParser & parser);
};


// forward declartions
struct DemultiplexingParams;
struct QualityTrimmingParams;
struct ProgramParams;
struct InputFileStreams;
struct ProcessingParams;

// structs
struct ProgramParams
{
    unsigned int fileCount;
    bool showSpeed;
    unsigned int firstReads;
    unsigned records;
    unsigned int num_threads;
    bool ordered;

    ProgramParams() : fileCount(0), showSpeed(false), firstReads(0), records(0), num_threads(0), ordered(false) {};
};

//Function declarations
seqan::ArgumentParser initParser(const FlexiProgram flexiProgram);

int loadBarcodes(char const * path, DemultiplexingParams& params);

int loadDemultiplexingParams(seqan::ArgumentParser const& parser, DemultiplexingParams& params);

int loadAdapterTrimmingParams(seqan::ArgumentParser const& parser, AdapterTrimmingParams & params);

int loadQualityTrimmingParams(seqan::ArgumentParser const & parser, QualityTrimmingParams & params);

int openStream(seqan::CharString const & file, seqan::SeqFileIn & inFile);

int loadProgramParams(seqan::ArgumentParser const & parser, ProgramParams& params, InputFileStreams& vars);

int checkParams(ProgramParams const & programParams, InputFileStreams const& inputFileStreams, ProcessingParams const & processingParams,
    DemultiplexingParams const & demultiplexingParams, AdapterTrimmingParams const & adapterTrimmingParams,
    QualityTrimmingParams & qualityTrimmingParams);

// function definitions
void ArgumentParserBuilder::addGeneralOptions(seqan::ArgumentParser & parser, const FlexiProgram flexiProgram)
{
    // GENERAL OPTIONS -----------------------------------------
    addSection(parser, "General Options");

    seqan::ArgParseOption showSpeedOpt = seqan::ArgParseOption(
        "ss", "showSpeed", "Show speed in base pairs per second.");
    addOption(parser, showSpeedOpt);

    seqan::ArgParseOption writeStatsOpt = seqan::ArgParseOption(
        "st", "writeStats", "Write statistics into a file");
    addOption(parser, writeStatsOpt);

    seqan::ArgParseOption recordOpt = seqan::ArgParseOption(
        "r", "records", "Number of records to be read in one run. Adjust this so that one batch of read can fit into your CPU cache.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(recordOpt, 1000);
    setMinValue(recordOpt, "1");
    addOption(parser, recordOpt);

    seqan::ArgParseOption noQualOpt = seqan::ArgParseOption(
        "nq", "noQualities", "Force .fa format for output files.");
    addOption(parser, noQualOpt);

    seqan::ArgParseOption orderedOpt = seqan::ArgParseOption(
        "od", "ordered", "Keep reads in order. Needs -r 1 option to work properly.");
    addOption(parser, orderedOpt);

    seqan::ArgParseOption firstReadsOpt = seqan::ArgParseOption(
        "fr", "reads", "Process only first n reads.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setMinValue(firstReadsOpt, "1");
    addOption(parser, firstReadsOpt);

    seqan::ArgParseOption threadOpt = seqan::ArgParseOption(
        "tnum", "threads", "Number of threads used.",
        seqan::ArgParseOption::INTEGER, "THREADS");
    setDefaultValue(threadOpt, 1);
    setMinValue(threadOpt, "1");
    addOption(parser, threadOpt);

    if (flexiProgram == FlexiProgram::ADAPTER_REMOVAL || flexiProgram == FlexiProgram::FILTERING || flexiProgram == FlexiProgram::QUALITY_CONTROL)
    {
        seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
            "o", "output", "Name of the output file.",
            seqan::ArgParseOption::OUTPUT_FILE, "OUTPUT");
        setValidValues(outputOpt, seqan::SeqFileOut::getFileExtensions());
        addOption(parser, outputOpt);
    }
    else
    {
        seqan::ArgParseOption outputOpt = seqan::ArgParseOption(
            "o", "output", "Prefix and file ending of output files (prefix$.fa - $: placeholder which will be determined by the program.).",
            seqan::ArgParseOption::OUTPUT_PREFIX, "OUTPUT");
        setValidValues(outputOpt, seqan::SeqFileOut::getFileExtensions());
        addOption(parser, outputOpt);
    }

    if (flexiProgram == FlexiProgram::ALL_STEPS)
    {
        seqan::ArgParseOption adTagOpt = seqan::ArgParseOption(
            "t", "tag", "Tags IDs of sequences which had adapters removed and/or were quality-trimmed.");
        addOption(parser, adTagOpt);
    }

    seqan::ArgParseOption adInfoOpt = seqan::ArgParseOption(
        "ni", "noInfo", "Don't print paramter overview to console.");
    addOption(parser, adInfoOpt);

    seqan::ArgParseOption finMinLenOpt = seqan::ArgParseOption(
        "fm", "finalMinLength", "Deletes read (and mate)"
        " if on of them is shorter than the given value after the complete worflow.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setMinValue(finMinLenOpt, "1");
    addOption(parser, finMinLenOpt);

    seqan::ArgParseOption finLenOpt = seqan::ArgParseOption(
        "fl", "finalLength", "Trims reads to desired length after the complete workflow.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(finLenOpt, 0);
    setMinValue(finLenOpt, "1");
    addOption(parser, finLenOpt);
}

void ArgumentParserBuilder::addFilteringOptions(seqan::ArgumentParser & parser)
{
    // Preprocessing and Filtering
    addSection(parser, "Preprocessing and Filtering");

    seqan::ArgParseOption leftTrimOpt = seqan::ArgParseOption(
        "tl", "trimLeft", "Number of Bases to be trimmed from the 5'end(s) before further processing.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(leftTrimOpt, 0);
    setMinValue(leftTrimOpt, "0");
    addOption(parser, leftTrimOpt);

    seqan::ArgParseOption tagTrimmingOpt = seqan::ArgParseOption(
        "tt", "tagTrimming", "Write trimmed-out segments into id.");
    addOption(parser, tagTrimmingOpt);

    seqan::ArgParseOption rigthTrimOpt = seqan::ArgParseOption(
        "tr", "trimRight", "Number of Bases to be trimmed from the 3'end(s) before further processing.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(rigthTrimOpt, 0);
    setMinValue(rigthTrimOpt, "0");
    addOption(parser, rigthTrimOpt);

    seqan::ArgParseOption minLenOpt = seqan::ArgParseOption(
        "ml", "minLength", "Required minimal length of reads after trimming.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(minLenOpt, 0);
    setMinValue(minLenOpt, "0");
    addOption(parser, minLenOpt);

    seqan::ArgParseOption uncalledOpt = seqan::ArgParseOption(
        "u", "uncalled", "Number of allowed uncalled bases per sequence.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(uncalledOpt, 0);
    setMinValue(uncalledOpt, "0");
    addOption(parser, uncalledOpt);

    seqan::ArgParseOption substituteOption = seqan::ArgParseOption(
        "s", "substitute", "Substitue Dna-character for N's.",
        seqan::ArgParseOption::STRING, "BASE");
    setDefaultValue(substituteOption, "A");
    setValidValues(substituteOption, "A C G T");
    addOption(parser, substituteOption);

    addTextSection(parser, "EXAMPLE");
    std::string appName;
    assign(appName, getAppName(parser));
    addText(parser, appName + " READS.fq -tr r -u 1 -o RESULT.fq");

}

void ArgumentParserBuilder::addDemultiplexingOptions(seqan::ArgumentParser & parser)
{
    // Barcode Demultiplexing
    addSection(parser, "Demultiplexing Options");

    seqan::ArgParseOption barcodeFileOpt = seqan::ArgParseOption(
        "b", "barcodes", "FastA file containing the used barcodes and their IDs. Necessary for demutiplexing.",
        seqan::ArgParseArgument::INPUT_FILE, "BARCODE_FILE");
    setValidValues(barcodeFileOpt, seqan::SeqFileIn::getFileExtensions());
    addOption(parser, barcodeFileOpt);

    seqan::ArgParseOption multiplexFileOpt = seqan::ArgParseOption(
        "x", "multiplex", "FastA/FastQ file containing the barcode for each read.",
        seqan::ArgParseArgument::INPUT_FILE, "MULTIPLEX_FILE");
    setValidValues(multiplexFileOpt, seqan::SeqFileIn::getFileExtensions());
    addOption(parser, multiplexFileOpt);

    addOption(parser, seqan::ArgParseOption(
        "app", "approximate", "Select approximate barcode demultiplexing, allowing one mismatch."));

    addOption(parser, seqan::ArgParseOption(
        "hc", "hardClip", "Select hardClip option, clipping the first length(barcode) bases in any case."));

    addOption(parser, seqan::ArgParseOption(
        "ex", "exclude", "Exclude unidentified reads from further processing."));

    addTextSection(parser, "EXAMPLE");
    std::string appName;
    assign(appName, getAppName(parser));
    addText(parser, appName + " READS.fq -b BARCODES.fa -o RESULT.fq");
}

void ArgumentParserBuilder::addAdapterTrimmingOptions(seqan::ArgumentParser & parser, const FlexiProgram flexiProgram)
{
    // ADAPTER TRIMMING
    addSection(parser, "Adapter removal options");
    seqan::ArgParseOption adapterFileOpt = seqan::ArgParseOption(
        "a", "adapters", "FastA file containing the two adapter sequences. "
        "The adapters according to the layout: 5'-adapter1-read-adapter2-3'.",
        seqan::ArgParseArgument::INPUT_FILE, "ADAPTER_FILE");
    setValidValues(adapterFileOpt, seqan::SeqFileIn::getFileExtensions());
    addOption(parser, adapterFileOpt);

    seqan::ArgParseOption noAdapterOpt = seqan::ArgParseOption(
        "pa", "paired-adapterTrimming", "Trim adapters from paired-end reads without using reference adapters.");
    addOption(parser, noAdapterOpt);

    seqan::ArgParseOption errorOpt = seqan::ArgParseOption(
        "e", "errors", "Allowed errors in adapter detection.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    addOption(parser, errorOpt);

    seqan::ArgParseOption rateOpt = seqan::ArgParseOption(
        "er", "error rate", "Allowed errors per overlap in adapter detection.",
        seqan::ArgParseOption::DOUBLE, "VALUE");
    setDefaultValue(rateOpt, 0.2);
    addOption(parser, rateOpt);

    seqan::ArgParseOption overlapOpt = seqan::ArgParseOption(
        "ol", "overlap", "Minimum length of overlap for a significant adapter match.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(overlapOpt, 4);
    addOption(parser, overlapOpt);

    seqan::ArgParseOption overhangOpt = seqan::ArgParseOption(
        "oh", "overhang", "Number of bases that the adapter can stick over at the opposite end",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(overhangOpt, 0);
    addOption(parser, overhangOpt);

    seqan::ArgParseOption timesOpt = seqan::ArgParseOption(
        "times", "times", "Do at maximum N iterations of adapter filtering. Every iteration the best matching adapter is removed.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(timesOpt, 1);
    addOption(parser, timesOpt);

    if (flexiProgram != FlexiProgram::ALL_STEPS)
    {
        seqan::ArgParseOption adTagOpt = seqan::ArgParseOption(
            "t", "tag", "Tags IDs of sequences which had adapters removed.");
        addOption(parser, adTagOpt);
    }

    addTextSection(parser, "EXAMPLE");
    std::string appName;
    assign(appName, getAppName(parser));
    addText(parser, appName + " READS.fq -a ADAPTER.fa -o RESULT.fq");

}

void ArgumentParserBuilder::addReadTrimmingOptions(seqan::ArgumentParser & parser, const FlexiProgram flexiProgram)
{
    // READ TRIMMING
    addSection(parser, "Quality trimming options");
    seqan::ArgParseOption qualOpt = seqan::ArgParseOption(
        "q", "quality", "Quality threshold for read trimming.",
        seqan::ArgParseArgument::INTEGER, "PHRED");
    setMinValue(qualOpt, "0");
    setMaxValue(qualOpt, "40");
    addOption(parser, qualOpt);

    seqan::ArgParseOption lenOpt = seqan::ArgParseOption(
        "l", "length",
        "Minimum read length after trimming. "
        "Shorter reads will be substituted by a single N or removed if the paired read is too short as well.",
        seqan::ArgParseArgument::INTEGER, "LENGTH");
    setDefaultValue(lenOpt, 1);
    setMinValue(lenOpt, "1");
    addOption(parser, lenOpt);

    seqan::ArgParseOption trimOpt = seqan::ArgParseOption(
        "m", "method", "Method for trimming reads.",
        seqan::ArgParseArgument::STRING, "METHOD");
    setDefaultValue(trimOpt, "WIN");
    setValidValues(trimOpt, "WIN BWA TAIL");
    addOption(parser, trimOpt);

    if (flexiProgram != FlexiProgram::ALL_STEPS)
    {
        seqan::ArgParseOption adTagOpt = seqan::ArgParseOption(
            "t", "tag", "Tags IDs of sequences which were quality-trimmed.");
        addOption(parser, adTagOpt);
    }

    addTextSection(parser, "EXAMPLE");
    std::string appName;
    assign(appName, getAppName(parser));
    addText(parser, appName + " READS.fq -q 20 -l 80 BARCODES -app -o RESULT.fq");
}

void FilteringParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Filtering Toolkit of Flexcat.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of Flexcat (a new implementation and extension of"
        " the original flexcat[1]) and can be used to filter reads and apply "
        "sequence independent trimming options");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
        "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
        "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "READS", true);
    setValidValues(fileArg, seqan::SeqFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

void AdapterRemovalParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The Adapter Removal Toolkit of Flexcat.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of Flexcat (a new implementation and extension of"
        " the original flexbar[1]) and removes adapter sequences from reads.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
        "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
        "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
#ifdef SEQAN_DATE
    seqan::setDate(parser, SEQAN_DATE);
#endif

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "READS", true);
    setValidValues(fileArg, seqan::SeqFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

void DemultiplexingParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Demultiplexing Toolkit of Flexcat.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of Flexcat (a new implementation and extension of"
        " the original flexbar[1]) and can be used for demultiplexing of reads.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
        "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
        "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");


    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "READS", true);
    setValidValues(fileArg, seqan::SeqFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

void QualityControlParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The SeqAn Quality Control Toolkit of Flexcat.");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "This program is a sub-routine of Flexcat (a new implementation and extension of"
        " the original flexbar [1]) and can be used for quality controlling of reads.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
        "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
        "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "READS", true);
    setValidValues(fileArg, seqan::SeqFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

void AllStepsParserBuilder::addHeader(seqan::ArgumentParser & parser)
{
    setCategory(parser, "NGS Quality Control");
    setShortDescription(parser, "The Flexcat NGS-Data Processing Toolkit");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI<[READ_FILE2]>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "Flexcat is a toolkit for the processing of sequenced NGS reads. "
        "It is a reimplementation and extension of the original Flexbar implementation of Dodt [1]. It is "
        "possible to demultiplex the reads and order them according to "
        "different kind of barcodes, to remove adapter contamination from "
        "reads, to trim low quality bases, filter N's or trim whole reads. The "
        "different tools are controlled through command line parameters and can "
        "operate on both single- and paired-end read data.");

    addDescription(parser, "[1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, "
        "C.  FLEXBAR—Flexible Barcode and Adapter Processing for "
        "Next-Generation Sequencing Platforms. Biology 2012, 1, 895-905.");

    addDescription(parser, "(c) Copyright 2015 by Benjamin Menkuec.");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "READS", true);
    setValidValues(fileArg, seqan::SeqFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "Either one (single-end) or two (paired-end) read files.");
}

// --------------------------------------------------------------------------
// Function initParser()
// --------------------------------------------------------------------------

//Defining the the argument parser
seqan::ArgumentParser initParser(const FlexiProgram flexiProgram)
{
    std::unique_ptr<ArgumentParserBuilder> argParseBuilder;

    switch (flexiProgram)
    {
    case FlexiProgram::ADAPTER_REMOVAL:
        argParseBuilder.reset(new AdapterRemovalParserBuilder());
        break;
    case FlexiProgram::DEMULTIPLEXING:
        argParseBuilder.reset(new DemultiplexingParserBuilder());
        break;
    case FlexiProgram::FILTERING:
        argParseBuilder.reset(new FilteringParserBuilder());
        break;
    case FlexiProgram::QUALITY_CONTROL:
        argParseBuilder.reset(new QualityControlParserBuilder());
        break;
    case FlexiProgram::ALL_STEPS:
        argParseBuilder.reset(new AllStepsParserBuilder());
        break;
    default:
        SEQAN_FAIL("Invalid program type.");
    }

    return argParseBuilder->build();
}

int loadBarcodes(char const * path, DemultiplexingParams& params)
{
    seqan::SeqFileIn bcFile;

    if (!open(bcFile, path, seqan::OPEN_RDONLY))
    {
        std::cerr << "Error while opening file'" << params.barcodeFile << "'.\n";
        return 1;
    }
    while (!atEnd(bcFile))
    {
        std::string id;
        std::string barcode;
        readRecord(id, barcode, bcFile);
        params.barcodeIds.emplace_back(id);
        params.barcodes.emplace_back(barcode);
    }
    if (params.approximate)                            //modifies the barcodes for approximate matching
    {
        buildAllVariations(params.barcodes);
    }
    return 0;
}

int loadDemultiplexingParams(seqan::ArgumentParser const& parser, DemultiplexingParams& params)
{
    // APPROXIMATE/EXACT MATCHING---------------------    
    params.approximate = seqan::isSet(parser, "app");
    // HARD CLIP MODE --------------------------------
    params.hardClip = seqan::isSet(parser, "hc") && !(isSet(parser, "ex"));
    // EXCLUDE UNIDENTIFIED --------------------------
    params.exclude = isSet(parser, "ex") && isSet(parser, "b");
    // RUN -------------------------------------------
    params.run = isSet(parser, "b");
    if (isSet(parser, "x"))
        params.runx = true;
    else
        params.runx = false;
    // BARCODES --------------------------------------
    if (isSet(parser, "b"))
    {
        getOptionValue(params.barcodeFile, parser, "b");
        if (loadBarcodes(seqan::toCString(params.barcodeFile), params) != 0)
            return 1;
    }
    return 0;
}

int loadAdapterTrimmingParams(seqan::ArgumentParser const& parser, AdapterTrimmingParams & params)
{
    // PAIRED-END ------------------------------
    int fileCount = getArgumentValueCount(parser, 0);
    // Only consider paired-end mode if two files are given and user wants paired mode.
    params.pairedNoAdapterFile = fileCount == 2 && isSet(parser, "pa");
    // Set run flag, depending on essential parameters.
    params.run = isSet(parser, "a") || (isSet(parser, "pa") && fileCount == 2);

    // Settings for adapter trimming
    int o;
    int e;
    double er;
    int oh;
    unsigned int times;
    getOptionValue(o, parser, "overlap");
    getOptionValue(e, parser, "e");
    getOptionValue(er, parser, "er");
    getOptionValue(oh, parser, "oh");
    getOptionValue(times, parser, "times");
    params.mode = AdapterMatchSettings(o, e, er, oh, times);

    // ADAPTER SEQUENCES ----------------------------
    std::string adapterFile, id;
    // If adapter sequences are given, we read them in any case.
    if (isSet(parser, "a"))
    {
        getOptionValue(adapterFile, parser, "a");
        seqan::SeqFileIn adapterInFile;
        if (!open(adapterInFile, seqan::toCString(adapterFile), seqan::OPEN_RDONLY))
        {
            std::cerr << "Error while opening file'" << adapterFile << "'.\n";
            return 1;
        }
        TAdapterSequence tempAdapter;
        AdapterItem adapterItem;
        unsigned int adapterId = 0;
        while (!atEnd(adapterInFile))
        {
            readRecord(id, adapterItem.seq, adapterInFile);
            if (id.find("3'") != std::string::npos)
                adapterItem.adapterEnd = AdapterItem::end3;
            else if (id.find("5'") != std::string::npos)
                adapterItem.adapterEnd = AdapterItem::end5;
            else
            {
                std::cerr << "End for adapter \"" << id << "\" not specified, adapter will be ignored.\n";
                continue;
            }
            adapterItem.anchored = id.find(":anchored:") != std::string::npos ? true : false;
            adapterItem.reverse = id.find(":reverse:") != std::string::npos ? true : false;

            adapterItem.overhang = oh;
            adapterItem.id = adapterId++;
            seqan::appendValue(params.adapters, adapterItem);
        }
    }
    // If they are not given, but we would need them (single-end trimming), output error.
    else if ((isSet(parser, "pa") && fileCount == 1))
    {
        std::cerr << "Unpaired adapter removal requires paired reads.\n";
        return 1;
    }
    else if (isSet(parser, "pa"))
    {
        std::cerr << "You can not specify adapters for paired adapter removal.\n";
        return 1;
    }

    return 0;
}

int loadQualityTrimmingParams(seqan::ArgumentParser const & parser, QualityTrimmingParams & params)
{
    // TRIMMING METHOD ----------------------------
    std::string method;
    getOptionValue(method, parser, "m");
    if (method == "WIN")
    {
        params.trim_mode = TrimmingMode::E_WINDOW;
    }
    else if (method == "BWA")
    {
        params.trim_mode = TrimmingMode::E_BWA;
    }
    else
    {
        params.trim_mode = TrimmingMode::E_TAIL;
    }
    // QUALITY CUTOFF ----------------------------
    if (isSet(parser, "q"))
    {
        getOptionValue(params.cutoff, parser, "q");
    }
    // MINIMUM SEQUENCE LENGTH -------------------
    getOptionValue(params.min_length, parser, "l");
    // Set run flag, depending on essential parameters (which are in a valid state at this point).
    params.run = isSet(parser, "q");
    return 0;
}

int openStream(seqan::CharString const & file, seqan::SeqFileIn & inFile)
{
    if (!open(inFile, seqan::toCString(file)))
    {
        std::cerr << "Error while opening input file '" << file << "'.\n";
        return 1;
    }
    return 0;
}

int loadProgramParams(seqan::ArgumentParser const & parser, ProgramParams& params, InputFileStreams& vars)
{
    params.fileCount = getArgumentValueCount(parser, 0);
    // Load files.
    seqan::CharString fileName1, fileName2;
    getArgumentValue(fileName1, parser, 0, 0);
    if (openStream(fileName1, vars.fileStream1) != 0)
    {
        return 1;
    }
    if (params.fileCount == 2)
    {
        getArgumentValue(fileName2, parser, 0, 1);
        if (openStream(fileName2, vars.fileStream2) != 0)
        {
            return 1;
        }
        if (value(format(vars.fileStream1)) != value(format(vars.fileStream2)))
        {
            std::cerr << "Input files must have the same file format.\n";
            return 1;
        }
    }
    params.showSpeed = isSet(parser, "ss");

    params.firstReads = std::numeric_limits<unsigned>::max();
    if (seqan::isSet(parser, "fr"))
        getOptionValue(params.firstReads, parser, "fr");

    params.num_threads = 1;
    getOptionValue(params.num_threads, parser, "tnum");
    omp_set_num_threads(params.num_threads);

    getOptionValue(params.records, parser, "r");
    getOptionValue(params.ordered, parser, "od");
    return 0;
}

int checkParams(ProgramParams const & programParams, InputFileStreams const& inputFileStreams, ProcessingParams const & processingParams,
    DemultiplexingParams const & demultiplexingParams, AdapterTrimmingParams const & adapterTrimmingParams,
    QualityTrimmingParams & qualityTrimmingParams)
{
    // Were there options that activated at least one processing stage?

    if (!(adapterTrimmingParams.run || qualityTrimmingParams.run || demultiplexingParams.run || processingParams.runPre
        || processingParams.runPost))
    {
        std::cerr << "\nNo processing stage was specified.\n";
        return 1;
    }
    // If quality trimming was selected, check if file format includes qualities.
    if (qualityTrimmingParams.run)
    {
        if ((value(format(inputFileStreams.fileStream1)) != seqan::Find<seqan::FileFormat<seqan::SeqFileIn>::Type, seqan::Fastq>::VALUE) ||
            ((programParams.fileCount == 2) &&
                (value(format(inputFileStreams.fileStream2)) != seqan::Find<seqan::FileFormat<seqan::SeqFileIn>::Type, seqan::Fastq>::VALUE)))
        {
            std::cerr << "\nQuality trimming requires quality information, please specify fq files." << std::endl;
            return 1;
        }
    }
    // If we don't demultiplex (and therefore take file names from the barcode IDs), set file names.
    if (!demultiplexingParams.run)
    {
        //std::cout << basename(toCString(fileName1)) << "\n\n";
    }
    return 0;
}




