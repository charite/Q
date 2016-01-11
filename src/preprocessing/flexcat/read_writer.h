// ==========================================================================
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================

#pragma once

#include <string>

class OutputStreams
{
    using TSeqStream = std::unique_ptr<seqan::SeqFileOut>;
    struct StreamPair
    {
        TSeqStream first;
        TSeqStream second;
        std::string firstFilename;
        std::string secondFilename;
    };
    //using TStreamPair = std::pair<TSeqStream, TSeqStream>;
    using TStreamPair = StreamPair;

    std::map<int, TStreamPair> fileStreams;
    const std::string basePath;
    std::string extension;

    template < typename TStream, template<typename> class TRead, typename TSeq,
        typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value  > >
        inline void writeRecord(TStream& stream, TRead<TSeq>&& read, bool = false)
    {
        seqan::writeRecord(*(stream.first), std::move(read.id), std::move(read.seq));
    }

    template <typename TStream, template<typename> class TRead, typename TSeq,
        typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value  > >
        inline void writeRecord(TStream& stream, TRead<TSeq>&& read)
    {
        seqan::writeRecord(*(stream.first), std::move(read.id), std::move(read.seq));
        seqan::writeRecord(*(stream.second), std::move(read.idRev), std::move(read.seqRev));
    }

    //Adds a new output streams to the collection of streams.
    std::string createStream(TSeqStream& stream, const std::string fileName, bool useDefault)
    {
        std::string path = getBaseFilename();
        if (fileName != "")
            path += "_";
        if (useDefault)
            path += "_result";

        path += fileName + extension;
        stream = std::make_unique<seqan::SeqFileOut>(path.c_str());
        return path;
    }


public:
    // The correct file extension is determined from the base path, according to the available
    // file extensions of the SeqFileOut and used for all stored files.
    OutputStreams(const std::string& base, bool /*noQuality*/) : basePath(base)
    {
        std::vector<std::string> tmpExtensions = seqan::SeqFileOut::getFileExtensions();
        for(const auto& tmpExtension : tmpExtensions)
        {
            if (seqan::endsWith(basePath, tmpExtension))
            {
                extension = tmpExtension;
                break;
            }
        }
    }

    inline std::string getBaseFilename(void) const
    {
        return prefix(basePath, length(basePath) - length(extension));
    }

    inline std::string getFilename(const int streamIndex) const
    {
        auto it = fileStreams.find(streamIndex);
        if (it != fileStreams.end())
            return it->second.firstFilename;
        return std::string();
    }

    void addStream(const std::string fileName, const int streamIndex, const bool useDefault)
    {
        auto filename = createStream(fileStreams[streamIndex].first, fileName, useDefault);
        fileStreams[streamIndex].firstFilename = filename;
    }
    
    void addStreams(const std::string fileName1, const std::string fileName2, const int streamIndex, const bool useDefault)
    {
        auto filename = createStream(fileStreams[streamIndex].first, fileName1, useDefault);
        fileStreams[streamIndex].firstFilename = filename;
        filename = createStream(fileStreams[streamIndex].second, fileName2, useDefault);
        fileStreams[streamIndex].secondFilename = filename;
    }

    void addStream(const std::string fileName, const int id)
    {
        addStream(fileName, id, false);
    }

    void addStreams(const std::string fileName1, const std::string fileName2, const int id)
    {
        addStreams(fileName1, fileName2, id, false);
    }

    //This method takes a String of integers and checks if these integers are
    //already associated with a stream. If not, a new stream is added and the opened
    //file is named according to the list of names. One or two files are created.
    template <typename TNames>
    void updateStreams(const TNames& names, const bool pair)
    {

        for (unsigned i = 0; i < length(names) + 1; ++i)
        {
            const unsigned streamIndex = i;
            // If no stream for this id exists, create one.
            if (fileStreams.find(streamIndex) == fileStreams.end())
            {
                // If the index is 0 (unidentified) create special stream.
                // Otherwise use index to get appropriate name for output file.
                std::string file;
                if (streamIndex > 0)
                    file = names[streamIndex - 1];
                else
                    file = "unidentified";
                // Add file extension to stream and create it.
                if (pair)
                {
                    std::string file2 = file;
                    file += "_1";
                    file2 += "_2";
                    addStreams(file, file2, streamIndex);
                }
                else
                {
                    addStream(file, streamIndex);
                }
            }
        }
    }

    template <template<typename> class TRead, typename TSeq, typename TNames>
    void writeSeqs(std::vector<TRead<TSeq>>&& reads, const TNames& names)
    {
        updateStreams(names, std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value);
        for(auto& read : reads)
        {
            const unsigned streamIndex = read.demuxResult;
            writeRecord(fileStreams[streamIndex], std::move(read));
        }
    }

    ~OutputStreams(){}

};


template<typename TOutputStreams, typename TProgramParams>
struct ReadWriter
{
private:
    //using TItem = std::tuple < std::unique_ptr<std::vector<TRead<TSeq>>>, decltype(DemultiplexingParams::barcodeIds), GeneralStats>;

    TOutputStreams& _outputStreams;
    const TProgramParams& _programParams;
    std::chrono::time_point<std::chrono::steady_clock> _startTime;
    std::chrono::time_point<std::chrono::steady_clock> _lastScreenUpdate;
    GeneralStats _stats;
public:
    ReadWriter(TOutputStreams& outputStreams, const TProgramParams& programParams) :
        _outputStreams(outputStreams), _programParams(programParams), _startTime(std::chrono::steady_clock::now()) {};

    template <typename TItem>
    void operator()(TItem item)
    {
        const auto t1 = std::chrono::steady_clock::now();
        _outputStreams.writeSeqs(std::move(*std::get<0>(*item)), std::get<1>(*item));
        _stats += std::get<2>(*item);

        // terminal output
        const auto ioTime = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - t1).count();
        _stats.ioTime += ioTime;
        const auto deltaLastScreenUpdate = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - _lastScreenUpdate).count();
        if (deltaLastScreenUpdate > 1)
        {
            const auto deltaTime = std::chrono::duration_cast<std::chrono::duration<float>>(std::chrono::steady_clock::now() - _startTime).count();
            if (_programParams.showSpeed)
                std::cout << "\rReads processed: " << _stats.readCount << "   (" << static_cast<int>(_stats.readCount / deltaTime) << " Reads/s)";
            else
                std::cout << "\rReads processed: " << _stats.readCount;
            _lastScreenUpdate = std::chrono::steady_clock::now();
        }
    }
    GeneralStats get_result()
    {
        return _stats;
    }
    void getStats(GeneralStats& stats)
    {
        stats = _stats;
    }
};
