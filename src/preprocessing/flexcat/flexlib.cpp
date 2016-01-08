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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// ==========================================================================

#include <algorithm>

#define SEQAN_PROFILE
#ifdef SEQAN_ENABLE_DEBUG
#define DEBUG_MSG(str) do { std::cerr << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#include <iostream>
#include <future>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#ifdef _WIN32
#include <direct.h>
#endif

// Headers for creating subdirectories.
#include <errno.h>

#include "flexcat.h"
#include "flexlib.h"
#include "argument_parser.h"
#include "read_trimming.h"
#include "adapter_trimming.h"
#include "demultiplex.h"
#include "general_processing.h"
#include "helper_functions.h"
#include "read.h"
#include "read_writer.h"
#include "ptc.h"

using namespace seqan;


// ============================================================================
// Functions
// ============================================================================


template<template <typename> class TRead, typename TSeq, typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value>>
inline void loadMultiplex(std::vector<TRead<TSeq>>& reads, unsigned records, seqan::SeqFileIn& multiplexFile)
{
    (void)reads;
    (void)multiplexFile;
    (void)records;
}

template<template <typename> class TRead, typename TSeq, typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value>>
inline void loadMultiplex(std::vector<TRead<TSeq>>& reads, unsigned records, seqan::SeqFileIn& multiplexFile, bool = false)
{
    seqan::String<char> id;
    unsigned int i = 0;
    while (i < records && !atEnd(multiplexFile))
    {
        readRecord(id, reads[i].demultiplex, multiplexFile);
        ++i;
    }
}

// PROGRAM STAGES ---------------------
//Preprocessing Stage
template<typename TReadSet, typename TStats>
void preprocessingStage(const ProcessingParams& processingParams, TReadSet& readSet, TStats& generalStats)
{
    if (processingParams.runPre)
    {
        //Trimming and filtering
        if (processingParams.trimLeft + processingParams.trimRight + processingParams.minLen != 0)
            preTrim(readSet, processingParams.trimLeft, processingParams.trimRight,
                processingParams.minLen, processingParams.tagTrimming, generalStats);
       // Detecting uncalled Bases
        if (processingParams.runCheckUncalled)
        {
            if (processingParams.runSubstitute)
                processN(readSet, processingParams.uncalled, processingParams.substitute, generalStats);
            else
                processN(readSet, processingParams.uncalled, NoSubstitute(), generalStats);
        }
    }
}

// DEMULTIPLEXING
template <typename TRead, typename TFinder, typename TStats>
int demultiplexingStage(const DemultiplexingParams& params, std::vector<TRead>& reads, TFinder& esaFinder,
    TStats& generalStats)
{
    if (!params.run)
        return 0;
    if (!params.approximate)
    {
        demultiplex(reads, esaFinder, params.hardClip, generalStats, ExactBarcodeMatching(), params.exclude);
    }
    else
    {
        if (!check(reads, params.barcodes, generalStats))            // On Errors with barcodes return 1;
            return 1;
        demultiplex(reads, esaFinder, params.hardClip, generalStats, ApproximateBarcodeMatching(), params.exclude);
    }
    return 0;
}

// ADAPTER TRIMMING
template <typename TRead, typename TStats>
void adapterTrimmingStage(const AdapterTrimmingParams& params, std::vector<TRead>& reads, TStats& stats)
{
    if (!params.run)
        return;
    if(params.tag)
        stripAdapterBatch(reads, params.adapters, params.mode, params.pairedNoAdapterFile, stats.adapterTrimmingStats, TagAdapter<true>());
    else
        stripAdapterBatch(reads, params.adapters, params.mode, params.pairedNoAdapterFile, stats.adapterTrimmingStats, TagAdapter<false>());
}

// QUALITY TRIMMING
template <typename TRead, typename TStats>
void qualityTrimmingStage(const QualityTrimmingParams& params, std::vector<TRead>& reads, TStats& stats)
{
    if (params.run)
    {
        switch (params.trim_mode)
        {
        case TrimmingMode::E_WINDOW:
        {
            trimBatch(reads, params.cutoff, Mean(5), params.tag);
            break;
        }
        case TrimmingMode::E_BWA:
        {
            trimBatch(reads, params.cutoff, BWA(), params.tag);
        }
        case TrimmingMode::E_TAIL:
        {
            trimBatch(reads, params.cutoff, Tail(), params.tag);
        }
        }
    }
    stats.removedQuality += removeShortSeqs(reads, params.min_length);
}

//Postprocessing
template<typename TRead, typename TStats>
void postprocessingStage(const ProcessingParams& params, std::vector<TRead>& reads, TStats& stats)
{
    if (params.runPost)
    {
        if ((params.finalMinLength != 0) && (params.finalLength == 0))
            stats.removedShort += removeShortSeqs(reads, params.finalMinLength);
        else if (params.finalLength != 0)
            trimTo(reads, params.finalLength, stats);
    }
}

// END PROGRAM STAGES ---------------------
template <typename TOutStream, typename TStats>
void printStatistics(const ProgramParams& programParams, const TStats& generalStats, DemultiplexingParams& demultiplexParams,
                const AdapterTrimmingParams& adapterParams, const OutputStreams& outputStreams, const bool timing, TOutStream &outStream)
{
    bool paired = programParams.fileCount == 2;
    bool adapter = adapterParams.run;
    outStream << std::endl;
    outStream << "Read statistics\n";
    outStream << "===============\n";
    outStream << "Reads processed:\t" << generalStats.readCount;
    if (paired) 
    {
        outStream << " (2 * " << generalStats.readCount << " single reads)";
    }
    outStream << std::endl;
    double dropped = generalStats.removedQuality + generalStats.removedN
        + generalStats.removedShort + (demultiplexParams.exclude * generalStats.removedDemultiplex);
    outStream << "  Reads dropped:\t" << dropped << "\t(" << std::setprecision(3) 
        << dropped / double(generalStats.readCount) * 100 << "%)\n";
    if (dropped != 0.0)
    {
        if (generalStats.removedN != 0)
        {
            outStream << "    Due to N content:\t" << generalStats.removedN << "\t(" << std::setprecision(3)
                << double(generalStats.removedN) / dropped * 100 << "%)\n";
        }
        if (generalStats.removedQuality != 0)
        {
            outStream << "    Due to quality:\t" << generalStats.removedQuality << "\t("
                << std::setprecision(3) << double(generalStats.removedQuality) / dropped * 100 << "%)\n";
        }
        if (generalStats.removedShort != 0)
        {
            outStream << "    Due to shortness:\t" << generalStats.removedShort << "\t("
                << std::setprecision(3) << double(generalStats.removedShort) / dropped * 100 << "%)\n";
        }
    }
    if (generalStats.uncalledBases != 0)
        {
            outStream << "  Remaining uncalled (or masked) bases: " << generalStats.uncalledBases << "\n";
        }
    //Statistics for Demultiplexing
    if (demultiplexParams.run)    
    {
        outStream << "\nBarcode Demultiplexing statistics\n";
        outStream << "=================================\n";

        unsigned barcodesTotal = length(demultiplexParams.barcodes);
        if (demultiplexParams.approximate)
        {
            barcodesTotal = barcodesTotal/(length(demultiplexParams.barcodes[0])*5);
        }

        outStream << "Reads per barcode:\n";
        outStream << "Unidentified:\t" << generalStats.matchedBarcodeReads[0];
        if (generalStats.readCount != 0)
        {
            outStream  << "\t\t(" << std::setprecision(3) << (double)generalStats.matchedBarcodeReads[0] /
                ((double)generalStats.readCount) * 100 << "%)";
        }  
        outStream << "\n";
        for (unsigned i = 1; i <= barcodesTotal; ++i)
        {
            outStream << demultiplexParams.barcodeIds[i-1]<<":\t" << generalStats.matchedBarcodeReads[i];
            if (generalStats.readCount != 0)
            {
                outStream  << "\t\t(" << std::setprecision(3) << (double)generalStats.matchedBarcodeReads[i] /
                    ((double)generalStats.readCount) * 100 << "%)";
            }
            outStream << "\n";
        }
        outStream << std::endl;
    }
    outStream << "File statistics\n";
    outStream << "===============\n";
    // How many reads are left.f
    int survived = generalStats.readCount - generalStats.removedN
        - generalStats.removedQuality - generalStats.removedShort
        - demultiplexParams.exclude * generalStats.removedDemultiplex;
    // In percentage points.
    double surv_proc = (double)survived / (double)generalStats.readCount * 100;
// TODO: add support for multiple output files
    outStream << outputStreams.getFilename(0) + ":\n";
    outStream << "-------\n";
    outStream << "  Surviving Reads: " << survived << "/" << generalStats.readCount
              << " (" << std::setprecision(3) << surv_proc << "%)\n";
    outStream << std::endl;
    if (adapter)
    {
        outStream << "ADAPTERS" << std::endl;
        outStream << "========" << std::endl;
        unsigned int i = 0;
        for (const auto adapterItem : adapterParams.adapters)
        {
            if (adapterItem.adapterEnd == AdapterItem::end3)
                outStream << " 3'";
            else
                outStream << " 5'";
            outStream << "-adapter, ";
            outStream << generalStats.adapterTrimmingStats.numRemoved[adapterItem.id] << "x, " << adapterItem.seq << std::endl;
            ++i;
        }
        outStream << std::endl;
        const auto totalRemoved = std::accumulate(generalStats.adapterTrimmingStats.numRemoved.begin(), generalStats.adapterTrimmingStats.numRemoved.end(), 0);
        outStream << "removed adapters: " << totalRemoved << "\n";
        outStream << std::endl;
        if (totalRemoved != 0)
        {
            int mean = generalStats.adapterTrimmingStats.overlapSum / totalRemoved;
            outStream << "Adapter sizes:\n";
            outStream << "Min: " << generalStats.adapterTrimmingStats.minOverlap << ", Mean: " << mean
                << ", Max: " << generalStats.adapterTrimmingStats.maxOverlap << "\n\n";
        }
        outStream << "Number of removed adapters\nmismatches\t0\t1\t2\t3\t4\t5\t6\t7\t8\nlength" << std::endl;
        i = 0;
        for (const auto adaptersSizeX : generalStats.adapterTrimmingStats.removedLength)
        {
            outStream << ++i<<"\t";
            for (const auto adaptersSizeXMismatchesN : adaptersSizeX)
                outStream << "\t" << adaptersSizeXMismatchesN;
            outStream << std::endl;
        }
    }

    // Print processing and IO time. IO is (approx.) the whole loop without the processing part.
    if (timing)
    {
        outStream << std::endl;
        outStream << "Time statistics:\n";
        outStream << "==================\n";
        outStream << "Processing time: " << std::setw(5) << generalStats.processTime << " seconds.\n";
        outStream << "       I/O time: " << std::setw(5) << generalStats.ioTime << " seconds.\n";
        outStream << std::endl;
    }
}

template < template<typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value> >
unsigned int readReads(std::vector<TRead<TSeq>>& reads, const unsigned int records, InputFileStreams& inputFileStreams, bool = false)
{
    reads.resize(records);
    unsigned int i = 0;
    TSeq seq;
    while (i < records && !atEnd(inputFileStreams.fileStream1))
    {
        readRecord(reads[i].id, seq, inputFileStreams.fileStream1);
        reads[i].seq = seq;
        ++i;
    }
    reads.resize(i);
    return i;
}

template < template<typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value> >
unsigned int readReads(std::vector<TRead<TSeq>>& reads, const unsigned int records, InputFileStreams& inputFileStreams)
{
    reads.resize(records);
    unsigned int i = 0;
    while (i < records && !atEnd(inputFileStreams.fileStream1))
    {
        readRecord(reads[i].id, reads[i].seq, inputFileStreams.fileStream1);
        readRecord(reads[i].idRev, reads[i].seqRev, inputFileStreams.fileStream2);
        ++i;
    }
    reads.resize(i);
    return i;
}

// END FUNCTION DEFINITIONS ---------------------------------------------
template<template <typename> class TRead, typename TSeq, typename TEsaFinder, typename TStats>
int mainLoop(TRead<TSeq>, const ProgramParams& programParams, InputFileStreams& inputFileStreams, const DemultiplexingParams& demultiplexingParams, const ProcessingParams& processingParams, const AdapterTrimmingParams& adapterTrimmingParams,
    const QualityTrimmingParams& qualityTrimmingParams, TEsaFinder& esaFinder,
    OutputStreams& outputStreams, TStats& stats)
{
    using TReadWriter = ReadWriter<OutputStreams, ProgramParams>;
    TReadWriter readWriter(outputStreams, programParams);

    unsigned int numReads = 0;
    auto readReader = [&numReads, &programParams, &inputFileStreams]() {
        auto item = std::make_unique<std::vector<TRead<TSeq>>>();
        if (numReads > programParams.firstReads)    // maximum read number reached -> dont do further reads
        {
            item.release();   // return empty unique_ptr to signal end
            return std::move(item);
        }
        readReads(*item, programParams.records, inputFileStreams);
        loadMultiplex(*item, programParams.records, inputFileStreams.fileStreamMultiplex);
        numReads += item->size();
        if (item->empty())    // no more reads available
            item.release();   // return empty unique_ptr to signal eof
        return std::move(item);
    };

    auto transformer = [&](auto reads){
        GeneralStats generalStats(length(demultiplexingParams.barcodeIds) + 1, adapterTrimmingParams.adapters.size());
        generalStats.readCount = reads->size();
        preprocessingStage(processingParams, *reads, generalStats);
        if (demultiplexingStage(demultiplexingParams, *reads, esaFinder, generalStats) != 0)
            std::cerr << "DemultiplexingStage error" << std::endl;
        adapterTrimmingStage(adapterTrimmingParams, *reads, generalStats);
        qualityTrimmingStage(qualityTrimmingParams, *reads, generalStats);
        postprocessingStage(processingParams, *reads, generalStats);
        return std::make_unique<std::tuple<decltype(reads), decltype(demultiplexingParams.barcodeIds), decltype(generalStats) >>(std::make_tuple(std::move(reads), demultiplexingParams.barcodeIds, generalStats));
    };

    TStats generalStats(length(demultiplexingParams.barcodeIds) + 1, adapterTrimmingParams.adapters.size());

    if (programParams.num_threads > 1)
    {
        if (programParams.ordered)
        {
            auto ptc_unit = ptc::ordered_ptc(readReader, transformer, readWriter, programParams.num_threads);
            ptc_unit->start();
            auto f = ptc_unit->get_future();
            stats = f.get();
        }
        else
        {
            auto ptc_unit = ptc::unordered_ptc(readReader, transformer, readWriter, programParams.num_threads);
            ptc_unit->start();
            auto f = ptc_unit->get_future();
            stats = f.get();
        }
    }
    else
    {
        std::unique_ptr<std::vector<TRead<TSeq>>> readSet;
        const auto tMain = std::chrono::steady_clock::now();
        while (generalStats.readCount < programParams.firstReads)
        {
            readSet.reset(new std::vector<TRead<TSeq>>(programParams.records));
            auto t1 = std::chrono::steady_clock::now();
            const auto numReadsRead = readReads(*readSet, programParams.records, inputFileStreams);
            generalStats.ioTime += std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - t1).count();
            if (numReadsRead == 0)
                break;

            auto res = transformer(std::move(readSet));
            generalStats += std::get<2>(*res);

            t1 = std::chrono::steady_clock::now();
            outputStreams.writeSeqs(std::move(*(std::get<0>(*res))), demultiplexingParams.barcodeIds);
            generalStats.ioTime += std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - t1).count();

            // Print information
            const auto deltaTime = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - tMain).count();
            if (programParams.showSpeed)
                std::cout << "\rreads processed: " << generalStats.readCount << "   (" << static_cast<int>(generalStats.readCount / deltaTime) << " Reads/s)";
            else
                std::cout << "\rreads processed: " << generalStats.readCount;
        }
        stats = generalStats;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------
// Program entry point.

int flexcatMain(const FlexiProgram flexiProgram, int argc, char const ** argv)
{
    SEQAN_PROTIMESTART(loopTime);
    seqan::ArgumentParser parser = initParser(flexiProgram);

    // Additional checks
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if input was successfully parsed.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check if one or two input files (single or paired-end) were given.
    int fileCount = getArgumentValueCount(parser, 0);
    if (!(fileCount == 1 || fileCount == 2)){
        printShortHelp(parser);
        return 1;
    }

    //--------------------------------------------------
    // Parse general parameters.
    //--------------------------------------------------

    seqan::CharString output;
    getOptionValue(output, parser, "output");

    //--------------------------------------------------
    // Parse pre- and postprocessing parameters.
    //--------------------------------------------------
    ProcessingParams processingParams;

    if(flexiProgram == FlexiProgram::FILTERING || flexiProgram == FlexiProgram::ALL_STEPS)
    {
        seqan::CharString substituteString;
        getOptionValue(substituteString, parser, "s");
        processingParams.substitute = substituteString[0];
        processingParams.runCheckUncalled = seqan::isSet(parser, "u");
        processingParams.runSubstitute = seqan::isSet(parser, "s");
        getOptionValue(processingParams.uncalled, parser, "u");
        if (seqan::isSet(parser, "s") && (!seqan::isSet(parser, "u") || processingParams.uncalled == 0))
        {
            std::cout << "\nWarning: Substitute for uncalled bases is set, but will be ineffective unless number of allowed"
                << " uncalled bases is set to a value higher than 0 by using the paramter -u [VALUE]\n";
            return 1;
        }
        getOptionValue(processingParams.trimLeft, parser, "tl");
        processingParams.tagTrimming = seqan::isSet(parser, "tt");
		getOptionValue(processingParams.trimRight, parser, "tr");
        getOptionValue(processingParams.minLen, parser, "ml");
        if (seqan::isSet(parser, "fl"))
        {
            getOptionValue(processingParams.finalLength , parser, "fl");
        }
        if (seqan::isSet(parser, "fm"))
        {
            getOptionValue(processingParams.finalMinLength , parser, "fm");
        }
        processingParams.runPre = ((processingParams.minLen + processingParams.trimLeft + processingParams.trimRight != 0)
            || isSet(parser, "u"));
        processingParams.runPost = (processingParams.finalLength + processingParams.finalMinLength != 0);
        if(flexiProgram == FlexiProgram::FILTERING)
        {
            processingParams.runPre = true;
            processingParams.runPost = true;
        }

    }
    //--------------------------------------------------
    // Parse demultiplexing parameters.
    //--------------------------------------------------

    DemultiplexingParams demultiplexingParams;
    seqan::SeqFileIn multiplexInFile;    //Initialising the SequenceStream for the multiplex file

    if(flexiProgram == FlexiProgram::DEMULTIPLEXING || flexiProgram == FlexiProgram::ALL_STEPS)
    {
        getOptionValue(demultiplexingParams.multiplexFile, parser, "x");
        if (isSet(parser, "x"))
        {
            if (!open(multiplexInFile, seqan::toCString(demultiplexingParams.multiplexFile)))
            {
                std::cerr << "Could not open file " << demultiplexingParams.multiplexFile << " for reading!" << std::endl;
                return 1;
            }
        }

        if (loadDemultiplexingParams(parser, demultiplexingParams) != 0)
            return 1;

        if(flexiProgram == FlexiProgram::DEMULTIPLEXING)
            demultiplexingParams.run = true;
    }

    //--------------------------------------------------
    // Process Barcodes
    //--------------------------------------------------

    BarcodeMatcher esaFinder(demultiplexingParams.barcodes);

    if(flexiProgram == FlexiProgram::DEMULTIPLEXING && (!isSet(parser, "x") && !isSet(parser, "b")))
    {
        std::cerr << "No Barcodefile was provided." << std::endl;
        return 1;
    }

    //--------------------------------------------------
    // Parse quality trimming parameters.
    //--------------------------------------------------

    QualityTrimmingParams qualityTrimmingParams;
    if(flexiProgram == FlexiProgram::QUALITY_CONTROL || flexiProgram == FlexiProgram::ALL_STEPS)
    {
        if ( loadQualityTrimmingParams(parser, qualityTrimmingParams) != 0 )
            return 1;

        if(flexiProgram == FlexiProgram::QUALITY_CONTROL)
            qualityTrimmingParams.run = true;
    }

    //--------------------------------------------------
    // Parse adapter trimming parameters.
    //--------------------------------------------------

    AdapterTrimmingParams adapterTrimmingParams;

    if(flexiProgram == FlexiProgram::ADAPTER_REMOVAL || flexiProgram == FlexiProgram::QUALITY_CONTROL|| flexiProgram == FlexiProgram::ALL_STEPS)
    {
        if(flexiProgram == FlexiProgram::ADAPTER_REMOVAL || flexiProgram == FlexiProgram::ALL_STEPS)
        {
            if (loadAdapterTrimmingParams(parser, adapterTrimmingParams) != 0)
            {
                return 1;
            }
        }
        const bool tagOpt = seqan::isSet(parser, "t");
        if (tagOpt && !(qualityTrimmingParams.run || adapterTrimmingParams.run))
        {
            std::cout << "\nWarning: Tag option was selected without selecting adapter removal or quality-trimming stage.\n";
                return 1;
        }
        adapterTrimmingParams.tag = qualityTrimmingParams.tag = tagOpt;
        if(flexiProgram == FlexiProgram::ADAPTER_REMOVAL)
            adapterTrimmingParams.run = true;
    }

    //--------------------------------------------------
    // General program parameters and additional checks.
    //--------------------------------------------------

    ProgramParams programParams;
    InputFileStreams inputFileStreams;
    if (loadProgramParams(parser, programParams, inputFileStreams) != 0)
        return 1;

    if (checkParams(programParams, inputFileStreams, processingParams, demultiplexingParams, adapterTrimmingParams, qualityTrimmingParams) != 0)
        return 1;

    //--------------------------------------------------
    // Processing
    //--------------------------------------------------

    // Prepare output stream object and initial mapping from StringSets to files.
    bool noQuality = (value(format(inputFileStreams.fileStream1)) ==
                      Find<FileFormat<seqan::SeqFileIn>::Type, Fasta>::VALUE);
    if (isSet(parser, "nq"))
    {
        noQuality = true;
    }

    seqan::CharString filename1;
    getArgumentValue(filename1, parser, 0, 0);

    bool useDefault = false;
    if (output == "")
    {
        output = filename1;
        useDefault = true;
    }

    OutputStreams outputStreams(seqan::toCString(output), noQuality);

    // Output additional Information on selected stages:
    if (!isSet(parser, "ni"))
    {
        seqan::CharString filename2, out;
        getOptionValue(out, parser, "output");
        std::cout << "Overview:\n";
        std::cout << "=========\n";
        std::cout << "Application Type: " << sizeof(void*) * 8 << " bit" << std::endl;
        getArgumentValue(filename1, parser, 0, 0);
        std::cout << "Forward-read file: " << filename1 << "\n";
        if (programParams.fileCount == 2)
        {
            getArgumentValue(filename2, parser, 0, 1);
            std::cout << "Backward-read file: " << filename2 << "\n";
        }
        else
        {
            std::cout << "Backward-read file: NONE\n";
        }
        if (isSet(parser, "output"))
        {
            std::cout << "Output-path: " << out <<"\n";
        }
        else
        {
            std::cout << "Output-path: Working Directory\n";
        }
        std:: cout << "\n"; 
        std::cout << "General Options:\n";
        std::cout << "\tReads per block: " << programParams.records * ((programParams.fileCount == 2) + 1)<< "\n";
        /*
        if (isSet(parser, "c"))
        {
            std::cout << "\tCompress output: YES" << "\n";
        }tag
        else
        {
            std::cout << "\tCompress output: NO" << "\n";
        }
        */
		if (isSet(parser, "fr") && (value(format(inputFileStreams.fileStream1)) !=
			Find<FileFormat<seqan::SeqFileIn>::Type, Fasta>::VALUE))
		{
			std::cout << "\tForce no-quality output: YES\n";
		}
		if (isSet(parser, "nq") && (value(format(inputFileStreams.fileStream1)) !=
                                    Find<FileFormat<seqan::SeqFileIn>::Type, Fasta>::VALUE))
        {
            std::cout << "\tProcess only first n reads: " << programParams.firstReads << "\n";
        }
        else if (value(format(inputFileStreams.fileStream1)) != Find<FileFormat<seqan::SeqFileIn>::Type, Fasta>::VALUE)
        {
            std::cout << "\tForce no-quality output: NO\n";
        }
        std::cout << "\tNumber of threads: " << programParams.num_threads << "\n";
        if (programParams.ordered)
            std::cout << "\tOrder policy: ordered" << std::endl;
        else
            std::cout << "\tOrder policy: unordered" << std::endl;
        if(flexiProgram == FlexiProgram::ADAPTER_REMOVAL || flexiProgram == FlexiProgram::QUALITY_CONTROL|| flexiProgram == FlexiProgram::ALL_STEPS)
        {
            if (isSet(parser, "t"))
            {
                std::cout << "\tTag quality-trimmed or adapter-removed reads: YES\n";
            }
            else if (adapterTrimmingParams.run||qualityTrimmingParams.run)
            {
                std::cout << "\tTag quality-trimmed or adapter-removed reads: NO\n";
            }
            std::cout << "\n";
        }
    
        if(flexiProgram == FlexiProgram::FILTERING || flexiProgram == FlexiProgram::ALL_STEPS)
        {
            std::cout << "Pre-, Postprocessing and Filtering:\n";
            std::cout << "\tPre-trim 5'-end length: " << processingParams.trimLeft << "\n";
			std::cout << "\tPre-trim 3'-end length: " << processingParams.trimRight << "\n";
            std::cout << "\tExclude reads shorter than: " << processingParams.minLen << "\n";
            if (isSet(parser, "u"))
            {
                std::cout << "\tAllowed uncalled bases per sequence: " << processingParams.uncalled << "\n";
            }
            if (isSet(parser, "s"))
            {
                std::cout << "\tSubstitute for uncalled bases: " << processingParams.substitute << "\n";
            }
            if (isSet(parser, "fm") && (processingParams.finalLength == 0))
            {
                std::cout << "\tMinimum sequence length after COMPLETE workflow: " << processingParams.finalMinLength << "\n";
            }
            else if (processingParams.finalLength != 0)
            {
                    std::cout << "\tTrim sequences after COMPLETE workflow to length: " <<processingParams.finalLength<< "\n";
            }
            std::cout << "\n";   
        }
        if (demultiplexingParams.run)
        {
            std::cout << "Barcode Demultiplexing:\n";
            std::cout << "\tBarcode file: " << demultiplexingParams.barcodeFile << "\n";
            if (demultiplexingParams.runx)
            {
                std::cout << "\tMultiplex barcode file: " << demultiplexingParams.multiplexFile << "\n";
            }
            else
            {
                std::cout << "\tMultiplex barcodes file:  NO" << demultiplexingParams.multiplexFile << "\n";
            }
            if (demultiplexingParams.approximate)
            {
                std::cout << "\tApproximate matching: YES\n";
            }
            else
            {
                std::cout << "\tApproximate matching: NO\n";
            }
            if (demultiplexingParams.hardClip && !demultiplexingParams.runx)
            {
                std::cout << "\tHardClip mode: YES\n";
            }
            else
            {
                std::cout << "\tHardClip mode: NO\n";
            }
            if (demultiplexingParams.exclude)
            {
                std::cout << "\tExclude unidentified sequences: YES\n";
            }
            else
            {
                std::cout << "\tExclude unidentified sequences: NO\n";
            }
            std::cout << "\n";
        }
        if (adapterTrimmingParams.run)
        {
            std::cout << "Adapter Removal:\n";
            if (isSet(parser, "a"))
            {
                seqan::CharString adapter;
                getOptionValue(adapter, parser, "a");
                std::cout << "\tAdapter file: " << adapter << "\n";
            }
            else
            {
                std::cout << "\tAdapter file: NONE\n";
            }
            if (programParams.fileCount == 2)
            {
                if (adapterTrimmingParams.pairedNoAdapterFile)
                {
                    std::cout << "\tPaired end Mode: without adapter file\n";
                }
                else
                {
                    std::cout << "\tPaired end Mode: with adapter file\n";
                }
            }
            if (isSet(parser, "e") && !isSet(parser, "er"))
            {
                unsigned e, o;
                getOptionValue(e, parser, "e");
                getOptionValue(o, parser, "overlap");
				std::cout << "\tAllowed mismatches: " << e << "\n";
                std::cout << "\tMinimum overlap " << o << "\n";
            }
			else if (!isSet(parser, "e") && isSet(parser, "er"))
			{
				unsigned o;
				double er;
				getOptionValue(er, parser, "er");
				getOptionValue(o, parser, "overlap");
				std::cout << "\tAllowed error rate: " << er << "\n";
				std::cout << "\tMinimum overlap " << o << "\n";
			}
			else
			{
				std::cout << "\nWarning: errors and error rate can not be specified both at the same time.\n";
				return 1;
			}
            if (isSet(parser, "oh"))
            {
                unsigned overhang;
                getOptionValue(overhang, parser, "oh");
                std::cout << "\tOverhang " << overhang << "\n";
            }
            unsigned times;
            getOptionValue(times, parser, "times");
            std::cout << "\tMax adapter trimming iterations " << times << "\n";
			std::cout << "\n";
        }
        if (qualityTrimmingParams.run)
        {
            std::cout << "Quality Trimming:\n";
            std::cout << "\tMinumum PHRED-Quality: " << qualityTrimmingParams.cutoff << "\n";
            std::string method;
            getOptionValue(method, parser, "m");
            std::cout << "\tMethod: " << method << "\n";
            if (isSet(parser, "l"))
            {
                std::cout << "\tMinimum length after trimming: " << qualityTrimmingParams.min_length << "\n";
            }
            std::cout << "\n";
        }
    }
    // Start processing. Different functions are needed for one or two input files.
    std::cout << "\nProcessing reads...\n" << std::endl;
    GeneralStats generalStats(length(demultiplexingParams.barcodeIds) + 1, adapterTrimmingParams.adapters.size());

    if (fileCount == 1)
    {
        if (!demultiplexingParams.run)
            outputStreams.addStream("", 0, useDefault);
        if(demultiplexingParams.runx)
            mainLoop(ReadMultiplex<seqan::Dna5QString>(), programParams, inputFileStreams, demultiplexingParams, processingParams, adapterTrimmingParams, qualityTrimmingParams, esaFinder, outputStreams, generalStats);
        else
            mainLoop(Read<seqan::Dna5QString>(), programParams, inputFileStreams, demultiplexingParams, processingParams, adapterTrimmingParams, qualityTrimmingParams, esaFinder, outputStreams, generalStats);
    }
     else
     {
         if (!demultiplexingParams.run)
             outputStreams.addStreams("", "", 0, useDefault);
         if (demultiplexingParams.runx)
             mainLoop(ReadMultiplexPairedEnd<seqan::Dna5QString>(), programParams, inputFileStreams, demultiplexingParams, processingParams, adapterTrimmingParams, qualityTrimmingParams, esaFinder, outputStreams, generalStats);
         else
             mainLoop(ReadPairedEnd<seqan::Dna5QString>(), programParams, inputFileStreams, demultiplexingParams, processingParams, adapterTrimmingParams, qualityTrimmingParams, esaFinder, outputStreams, generalStats);
    }
    double loop = SEQAN_PROTIMEDIFF(loopTime);
    generalStats.processTime = loop - generalStats.ioTime;

    printStatistics(programParams, generalStats, demultiplexingParams, adapterTrimmingParams, outputStreams, !isSet(parser, "ni"), std::cout);
    if (isSet(parser, "st"))
    {
        std::fstream statFile;
#ifdef _MSC_VER
        statFile.open(std::string(seqan::toCString(outputStreams.getBaseFilename())) + "_flexcat_statistics.txt", std::fstream::out, _SH_DENYNO);
#else
        statFile.open(std::string(seqan::toCString(outputStreams.getBaseFilename())) + "_flexcat_statistics.txt", std::fstream::out);
#endif
        statFile << "command line: ";
        for (int i = 0;i < argc;++i)
            statFile << argv[i] << " ";
        statFile << std::endl;
        printStatistics(programParams, generalStats, demultiplexingParams, adapterTrimmingParams, outputStreams, !isSet(parser, "ni"), statFile);
        statFile.close();
    }
    return 0;
}
