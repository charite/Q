/*
Author: Benjamin Menkuec
Copyright 2015 Benjamin Menkuec
License: LGPL
*/

#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <string>
#include <algorithm>
#include <chrono>
#include <cassert>

#include "peak.h"
#include "BamRecordKey.h"

struct Statistics
{
    unsigned totalReads = 0;
    unsigned totalMappedReads = 0;
    unsigned filteredReads = 0;
    unsigned removedReads = 0;
    unsigned totalSamePositionReads = 0;
    unsigned readsAfterFiltering = 0;
    unsigned int couldNotMap = 0;
    unsigned int couldNotMapUniquely = 0;
};

template <typename TStream>
void printStatistics(TStream &stream, const Statistics &stats, const bool clusterFiltering, const bool tabbed=false)
{
    if (tabbed)
    {
        stream << "Total reads" << "\t"<< stats.totalReads << std::endl;
        stream << "Filtered reads (-fc Option)" << "\t" << stats.filteredReads << std::endl;
        stream << "Mapped reads" << "\t" << stats.totalMappedReads << std::endl;
        stream << "Non mappable reads" << "\t" << stats.couldNotMap << std::endl;
        stream << "Non uniquely mappable reads" << "\t" << stats.couldNotMapUniquely << std::endl;
        stream << "After barcode filtering" << "\t" << stats.totalMappedReads - stats.removedReads << "\t" << " (-" << stats.removedReads << ")" << std::endl;
        stream << "PCR duplication rate" << "\t" << static_cast<float>(stats.removedReads) / static_cast<float>(stats.totalMappedReads) << std::endl;
        stream << "Total duplet reads" << "\t" << stats.totalSamePositionReads << std::endl;
        if (clusterFiltering)
            stream << "After cluster filtering" << "\t" << stats.readsAfterFiltering << "\t" << " (-" << stats.totalMappedReads - stats.removedReads - stats.readsAfterFiltering << ")" << std::endl;
    }
    else
    {
        stream << "Total reads                          : " << stats.totalReads << std::endl;
        stream << "Filtered reads (-fc Option)          : " << stats.filteredReads << std::endl;
        stream << "Mapped reads                         : " << stats.totalMappedReads << std::endl;
        stream << "Non mappable reads                   : " << stats.couldNotMap << std::endl;
        stream << "Non uniquely mappable reads          : " << stats.couldNotMapUniquely << std::endl;
        stream << "After barcode filtering              : " << stats.totalMappedReads - stats.removedReads << " (-" << stats.removedReads << ")" << std::endl;
        stream << "PCR duplication rate                 : " << static_cast<float>(stats.removedReads) / static_cast<float>(stats.totalMappedReads) << std::endl;
        stream << "Total duplet reads                   : " << stats.totalSamePositionReads << std::endl;
        if (clusterFiltering)
            stream << "After cluster filtering              : " << stats.readsAfterFiltering << " (-" << stats.totalMappedReads - stats.removedReads - stats.readsAfterFiltering << ")" << std::endl;
    }
}

seqan::ArgumentParser buildParser(void)
{
    seqan::ArgumentParser parser;

    setCategory(parser, "Chip Nexus/Exo Data Processing");
    setShortDescription(parser, "Preprocessing Pipeline for Chip-Nexus and Chip-Exo data");
    addUsageLine(parser, " \\fI<READ_FILE1>\\fP \\fI[OPTIONS]\\fP");
    addDescription(parser,
        "");

    addDescription(parser, "");

    seqan::setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE, "mapped reads", true);
    setValidValues(fileArg, seqan::BamFileIn::getFileExtensions());
    addArgument(parser, fileArg);
    setHelpText(parser, 0, "SAM or BAM file");

    seqan::ArgParseOption recordFilterCluster = seqan::ArgParseOption(
        "f", "cluster size", "Minimum number of mapped reads at the same position",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(recordFilterCluster, 0);
    setMinValue(recordFilterCluster, "0");
    addOption(parser, recordFilterCluster);

    seqan::ArgParseOption recordWriteBed = seqan::ArgParseOption(
        "b", "bedGraph", "Create a BedGraph file");
    addOption(parser, recordWriteBed);

    seqan::ArgParseOption outputArtifactsOpt = seqan::ArgParseOption(
        "oa", "outputArtifacts", "Write PCR artifacts to BAM file");
    addOption(parser, outputArtifactsOpt);

    seqan::ArgParseOption randomSplitOpt = seqan::ArgParseOption(
        "rs", "randomSplit", "Split output randomly into two files");
    addOption(parser, randomSplitOpt);

    seqan::ArgParseOption recordOpt = seqan::ArgParseOption(
        "r", "records", "Number of records to be read in one run.",
        seqan::ArgParseOption::INTEGER, "VALUE");
    setDefaultValue(recordOpt, 10000);
    setMinValue(recordOpt, "10");
    addOption(parser, recordOpt);

    seqan::ArgParseOption filterChromosomesOpt = seqan::ArgParseOption(
        "fc", "filterChromosomes", "Regular expression to remove chromosomes",
        seqan::ArgParseOption::STRING, "REGEX");
    addOption(parser, filterChromosomesOpt);

    return parser;
}

std::string getFilePath(const std::string& fileName)
{
    std::size_t found = fileName.find_last_of("/\\");
    if (found == std::string::npos || found < 1)
        return std::string();
    return fileName.substr(0, found);
}

std::string getFilePrefix(const std::string& fileName, const bool withPath = true)
{
    std::size_t found = fileName.find_last_of(".");
    std::size_t found2 = std::string::npos;
    if (!withPath)
        found2 = fileName.find_last_of("/\\");
    if (found == std::string::npos)
        return std::string();
    if (found2 == std::string::npos)
        return fileName.substr(0, found);
    return fileName.substr(found2 + 1, found - found2);
}

template <typename TContext>
struct SaveBam
{
    SaveBam(const seqan::BamHeader header, TContext& context, const std::string& filename)
        : bamFileOut(static_cast<TContext>(context))
    {
        if (!open(bamFileOut, (filename + ".bam").c_str()))
        {
            std::cerr << "ERROR: Could not open " << filename << " for writing.\n";
            return;
        }
        writeHeader(bamFileOut, header);
    }
    void write(const seqan::BamAlignmentRecord& record)
    {
        writeRecord(bamFileOut, record);
    }
    void close()
    {
        seqan::close(bamFileOut);
    }
    seqan::BamFileOut bamFileOut;
};

typedef std::map<BamRecordKey<NoBarcode>, std::pair<unsigned int, unsigned int>> OccurenceMap;

BamRecordKey<NoBarcode> getKey(const OccurenceMap::value_type& val)
{
    return val.first;
}

unsigned getUniqueFrequency(const OccurenceMap::value_type& val)
{
    return val.second.second;
}

template <typename TOccurenceMap, typename TArtifactWriter, typename TChromosomeFilter, typename TBamWriter>
void processBamFile(seqan::BamFileIn& bamFileIn, const TArtifactWriter& artifactWriter, const TBamWriter& bamWriter, 
    const TChromosomeFilter& chromosomeFilter, TOccurenceMap &occurenceMap, Statistics& stats)
{
    seqan::BamAlignmentRecord record;
    unsigned tagID = 0;
    std::set<BamRecordKey<WithBarcode>> keySet;

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        if (atEnd(bamFileIn))
            break;
        ++stats.totalReads;
        if (chromosomeFilter.find(record.rID) != chromosomeFilter.end())
        {
            ++stats.filteredReads;
            continue;
        }
        const seqan::BamTagsDict tags(record.tags);
        if (seqan::findTagKey(tagID, tags, seqan::CharString("XM")))
        {
            __int32 tagValue = 0;
            extractTagValue(tagValue, tags, tagID);
            if (tagValue == 0)
                ++stats.couldNotMap;
            else
                ++stats.couldNotMapUniquely;
        }
        if (record.flag != 0x00 && record.flag != 0x10)
            continue;

        ++stats.totalMappedReads;
        const BamRecordKey<WithBarcode> key(record);
        const BamRecordKey<NoBarcode> pos(record);
        const auto insertResult = keySet.insert(key);
        OccurenceMap::mapped_type &mapItem = occurenceMap[pos];
        if (!insertResult.second)  // element was not inserted because it existed already
        {
            artifactWriter(std::move(record));
            ++stats.removedReads; // stats.removedReads = total non unique hits
        }
        else
        {
            bamWriter(std::move(record));
            ++mapItem.second; // unique hits
        }
        // total hits
        ++mapItem.first;
    }
}

int main(int argc, char const * argv[])
{
    // Additional checks
    seqan::ArgumentParser parser = buildParser();
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if input was successfully parsed.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check if one input file is given.
    int fileCount = getArgumentValueCount(parser, 0);
    if (fileCount != 1) {
        printShortHelp(parser);
        return 1;
    }

    unsigned numRecords;
    getOptionValue(numRecords, parser, "r");

    seqan::CharString fileName1, fileName2;
    getArgumentValue(fileName1, parser, 0, 0);

    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn(seqan::toCString(fileName1));
    
    std::string outFilename;
    outFilename = getFilePrefix(seqan::toCString(fileName1)) + std::string("_filtered");

    const bool filter = seqan::isSet(parser, "f");
    const bool bedOutputEnabled = seqan::isSet(parser, "b");
    const bool outputArtifacts = seqan::isSet(parser, "oa");
    const bool randomSplit = seqan::isSet(parser, "rs");
    seqan::CharString _filterChromosomes;
    seqan::getOptionValue(_filterChromosomes, parser, "fc");
    std::string filterChromosomes = seqan::toCString(_filterChromosomes);

    OccurenceMap occurenceMap;
    Statistics stats;

    std::cout << "barcode filtering... ";
    auto t1 = std::chrono::steady_clock::now();
    seqan::BamAlignmentRecord record;
    seqan::BamHeader header;
    readHeader(header, bamFileIn);
    const auto chromosomeFilterSet = calculateChromosomeFilter(filterChromosomes, contigNames(context(bamFileIn)));
    std::vector<seqan::BamAlignmentRecord> artifacts;
    srand(time(NULL));
    // dynamic parameter dispatching to processBamFile depending on randomSplit and outputArtifacts
    auto artifactWriter = [&artifacts](seqan::BamAlignmentRecord&& record) {return artifacts.emplace_back(record);};
    auto noArtifactWriter = [](seqan::BamAlignmentRecord&& record) {(void)record;return;};
    if (randomSplit)
    {
        SaveBam<seqan::BamFileIn> saveBam(header, bamFileIn, outFilename);
        SaveBam<seqan::BamFileIn> saveBamSplit1(header, bamFileIn, outFilename + "_split1");
        SaveBam<seqan::BamFileIn> saveBamSplit2(header, bamFileIn, outFilename + "_split2");
        auto bamWriterSplit = [&saveBam, &saveBamSplit1, &saveBamSplit2](seqan::BamAlignmentRecord&& record) {
            saveBam.write(record);
            if (rand() % 2)
                saveBamSplit1.write(record);
            else
                saveBamSplit2.write(record);
            return;};

        if(outputArtifacts)
            processBamFile(bamFileIn, artifactWriter, bamWriterSplit, chromosomeFilterSet, occurenceMap, stats);
        else
            processBamFile(bamFileIn, noArtifactWriter, bamWriterSplit, chromosomeFilterSet, occurenceMap, stats);
        saveBam.close();
        saveBamSplit1.close();
        saveBamSplit2.close();
    }
    else
    {
        SaveBam<seqan::BamFileIn> saveBam(header, bamFileIn, outFilename);
        auto bamWriter = [&saveBam](seqan::BamAlignmentRecord&& record) {return saveBam.write(record);};

        if (outputArtifacts)
            processBamFile(bamFileIn, artifactWriter, bamWriter, chromosomeFilterSet, occurenceMap, stats);
        else
            processBamFile(bamFileIn, noArtifactWriter, bamWriter, chromosomeFilterSet, occurenceMap, stats);
        saveBam.close();
    }
    auto t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;

    if (outputArtifacts)
    {
        std::cout << "writing artifacts to file... ";
        auto t1 = std::chrono::steady_clock::now();
        SaveBam<seqan::BamFileIn> saveArtifactsBam(header, bamFileIn, getFilePrefix(seqan::toCString(fileName1)) + "_artifacts");
        for (auto element : artifacts)
            saveArtifactsBam.write(std::move(element));
        saveArtifactsBam.close();
        auto t2 = std::chrono::steady_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    }

    const std::string outFilename2 = getFilePrefix(seqan::toCString(fileName1)) + std::string("_filtered2");
    if (filter)
    {
        t1 = std::chrono::steady_clock::now();
        std::cout << "filtering reads... ";
        unsigned clusterSize = 0;
        getOptionValue(clusterSize, parser, "f");
        seqan::BamFileIn bamFileIn2(seqan::toCString(outFilename + ".bam"));
        clear(header); // this is a bug in seqan, header has to be cleared before call to readHeader
        readHeader(header, bamFileIn2);
        SaveBam<seqan::BamFileIn> saveBam2(header, bamFileIn2, outFilename2);
        while (!atEnd(bamFileIn2))
        {
            readRecord(record, bamFileIn2);
            const auto occurenceMapIt =  occurenceMap.find(BamRecordKey<NoBarcode>(record));
            if (occurenceMapIt->second.second >= clusterSize)
            {
                saveBam2.write(record);
                ++stats.readsAfterFiltering;
            }
        }
        saveBam2.close();

        t2 = std::chrono::steady_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    }

    // occurenceMap = number of mappings for each location in genome
    // duplicationRate[x] = number of locations with x mappings
    t1 = std::chrono::steady_clock::now();
    std::cout << "calculating unique/non unique duplication Rate... ";

    SaveBed<seqan::BedRecord<seqan::Bed4>> saveBedForwardStrand(outFilename + "_forward");
    SaveBed<seqan::BedRecord<seqan::Bed4>> saveBedReverseStrand(outFilename + "_reverse");
    if (bedOutputEnabled)
    {
        saveBedForwardStrand.writeHeader("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
        saveBedReverseStrand.writeHeader("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
    }

    seqan::BedRecord<seqan::Bed4> bedRecord;
    std::vector<unsigned> duplicationRateUnique;
    std::vector<unsigned> duplicationRate;
    std::for_each(occurenceMap.begin(), occurenceMap.end(), [&](const OccurenceMap::value_type& val)
    {
        if (bedOutputEnabled)
        {
            bedRecord.rID = getKey(val).getRID();
            bedRecord.ref = contigNames(bamFileIn.context)[bedRecord.rID];
            if (getKey(val).isReverseStrand())    // reverse strand
            {
                // I think this -1 is not neccessary, but its here to reproduce the data from the CHipNexus paper exactly
                bedRecord.beginPos = getKey(val).get5EndPosition();
                bedRecord.endPos = bedRecord.beginPos + 1;
                bedRecord.name = std::to_string(-static_cast<int32_t>(val.second.second)); // abuse name as val parameter in BedGraph
                saveBedReverseStrand.write(bedRecord);
            }
            else    // forward strand
            {
                bedRecord.beginPos = getKey(val).get5EndPosition();
                bedRecord.endPos = bedRecord.beginPos + 1;
                bedRecord.name = std::to_string(val.second.second); // abuse name as val parameter in BedGraph
                saveBedForwardStrand.write(bedRecord);
            }
        }
        const auto uniqueHits = val.second.second;
        const auto totalHits = val.second.first;
        if (duplicationRateUnique.size() < uniqueHits)
            duplicationRateUnique.resize(uniqueHits);
        ++duplicationRateUnique[uniqueHits - 1];
        if (totalHits - uniqueHits > 0)
        {
            const auto nonUniqueHits = totalHits - uniqueHits;
            if (duplicationRate.size() < nonUniqueHits)
                duplicationRate.resize(nonUniqueHits);
            ++duplicationRate[nonUniqueHits - 1];
        }
    });
    saveBedForwardStrand.close();
    saveBedReverseStrand.close();

    std::fstream fs,fs2,fs3;
#ifdef _MSC_VER
    fs.open(getFilePrefix(seqan::toCString(fileName1)) + "_duplication_rate_positions.txt", std::fstream::out, _SH_DENYNO);
    fs2.open(getFilePrefix(seqan::toCString(fileName1)) + "_duplication_rate_reads.txt", std::fstream::out, _SH_DENYNO);
#else
    fs.open(getFilePrefix(seqan::toCString(fileName1)) + "_duplication_rate_positions.txt", std::fstream::out);
    fs2.open(getFilePrefix(seqan::toCString(fileName1)) + "_duplication_rate_reads.txt", std::fstream::out);
#endif
    fs << "rate" << "\t" << "unique" << "\t" << "non_unique" << std::endl;
    fs2 << "rate" << "\t" << "unique" << "\t" << "non_unique" << std::endl;
    unsigned maxLen = duplicationRateUnique.size() > duplicationRate.size() ? duplicationRateUnique.size() : duplicationRate.size();
    duplicationRateUnique.resize(maxLen);
    duplicationRate.resize(maxLen);
    std::vector<unsigned>::iterator it = duplicationRateUnique.begin();
    unsigned int sumUniqueReads = 0;
    unsigned int sumNonUniqueReads = 0;
    for (unsigned i = 0; i < maxLen;++i)
    {
        if (i > 0)
            stats.totalSamePositionReads += (duplicationRate[i] * (i+1));
        fs << i + 1 << "\t" << duplicationRateUnique[i] << "\t" << duplicationRate[i] << std::endl;
        fs2 << i + 1 << "\t" << duplicationRateUnique[i]*(i+1) << "\t" << (duplicationRate[i])*(i+1) <<std::endl;
        sumUniqueReads += duplicationRateUnique[i] * (i + 1);
        sumNonUniqueReads += duplicationRate[i] * (i+1);
    }
    assert(sumUniqueReads == stats.totalMappedReads - stats.removedReads);
    assert(stats.totalReads == stats.couldNotMap + stats.couldNotMapUniquely + stats.totalMappedReads + stats.filteredReads);
    assert(sumNonUniqueReads == stats.removedReads);
    t2 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << "s" << std::endl;
    duplicationRateUnique.clear();
    duplicationRate.clear();

    printStatistics(std::cout, stats, seqan::isSet(parser, "f"));
#ifdef _MSV_VER
    fs3.open(getFilePrefix(seqan::toCString(fileName1)) + "_nexcat_statistics.txt", std::fstream::out, _SH_DENYNO);
#else
    fs3.open(getFilePrefix(seqan::toCString(fileName1)) + "_nexcat_statistics.txt", std::fstream::out);
#endif
    printStatistics(fs3, stats, seqan::isSet(parser, "f"), true);
    fs3 << "Command line\t";
    for (int n = 0; n < argc;++n)
        fs3 << argv[n] << "\t";
    fs3.close();
	return 0;
}
