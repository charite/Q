// ==========================================================================
//                             adapterTrimming.h
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

#pragma once

#include <seqan/align.h>
#include "helper_functions.h"
#include "general_stats.h"

// ============================================================================
// Metafunctions
// ============================================================================


//brief A metafunction which constructs a ModifiedString type for an Alphabet.
template <class TValue>
struct STRING_REVERSE_COMPLEMENT
{
    typedef seqan::ModifiedString<
        seqan::ModifiedString<	seqan::String<TValue>, seqan::ModView<seqan::FunctorComplement<TValue> > >,
        seqan::ModReverse>	Type;
};


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Custom scoring matrix for ignoring alignments with N.
namespace seqan{
	//Struct used to define a new custom scoring matrix.
	struct AdapterScoringMatrix {};
     //Struct containing data for the seqan::AdapterScoringMatrix custom scoring matrix.
	 //Matches score 1, mismatches -1 and matches against N with 0.
	template <>
	struct ScoringMatrixData_<int, Dna5, AdapterScoringMatrix> {
		enum {
			VALUE_SIZE = ValueSize<Dna5>::VALUE,
			TAB_SIZE = VALUE_SIZE * VALUE_SIZE
		};
		static inline int const * getData() {
			static int const _data[TAB_SIZE] = {
			   1, -1, -1, -1, 0,
			  -1,  1, -1, -1, 0,
			  -1, -1,  1, -1, 0,
			  -1, -1, -1,  1, 0,
			   0,  0,  0,  0, 0
			};
			return _data;
		}
	};
}

using TAdapterAlphabet = seqan::Dna5Q;
using TAdapterSequence = seqan::String<TAdapterAlphabet>;

struct AdapterItem;
typedef std::vector< AdapterItem > AdapterSet;
using TReverseComplement = STRING_REVERSE_COMPLEMENT<TAdapterAlphabet>::Type;

struct AdapterItem
{
    using TAdapterSequence = seqan::String<TAdapterAlphabet>;

    enum AdapterEnd
    {
        end3,
        end5,
    };
    AdapterEnd adapterEnd;
    unsigned int overhang;
    unsigned id;
    bool anchored;
    bool reverse;
    TAdapterSequence seq;

    AdapterItem() : adapterEnd(end3), overhang(0), id(0), anchored(false), reverse(false){};
    AdapterItem(const TAdapterSequence &adapter) : adapterEnd(end3), overhang(0), id(0), anchored(false), reverse(false), seq(adapter){};
    AdapterItem(const TAdapterSequence &adapter, const AdapterEnd adapterEnd, const unsigned overhang, const unsigned id, const bool anchored, const bool reverse)
        : adapterEnd(adapterEnd), overhang(overhang), id(id), anchored(anchored), reverse(reverse), seq(adapter) {};

    AdapterItem getReverseComplement() const noexcept
    {
        auto seqCopy = seq;
        return AdapterItem(TReverseComplement(seqCopy), adapterEnd, overhang, id, anchored, reverse);
    }


    // the anchored matching mode should have the same functionality as cutadapt's anchored mode
    // see http://cutadapt.readthedocs.org/en/latest/guide.html
    // In anchored mode, the adapter has to start (3" adapter) or has to end (5" adapter) with the sequence
    // The anchored mode is rarely used, at least for 3" adapters
    // Todo: implement rooted mode

};

// Define scoring function type.
typedef seqan::Score<int, seqan::ScoreMatrix<seqan::Dna5, seqan::AdapterScoringMatrix> > TScore;

struct AdapterMatchSettings
{
    AdapterMatchSettings(const int m, const int e, const double er, const unsigned int oh, const unsigned int times) : min_length(m), errors(e), overhang(oh), errorRate(er), times(times)
    {}
    AdapterMatchSettings() : min_length(0), errors(0), overhang(0), errorRate(0), times(1) {};
   
    unsigned int min_length; //The minimum length of the overlap.
	int errors;     //The maximum number of errors we allow.
    unsigned int overhang;
	double errorRate;  //The maximum number of errors allowed per overlap
    unsigned int times;
};


struct AdapterTrimmingParams
{
    bool pairedNoAdapterFile;
    bool run;
    AdapterSet adapters;
    AdapterMatchSettings mode;
    bool tag;
    AdapterTrimmingParams() : pairedNoAdapterFile(false), run(false), tag(false) {};
};

// ============================================================================
// Functions
// ============================================================================


template <typename TRow>
unsigned countTotalGaps(TRow& row) noexcept
{
	return length(row) - length(source(row));
}

template <typename TAlign>
unsigned getOverlap(TAlign& align) noexcept
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	const TRow &row1 = seqan::row(align,0);
	const TRow &row2 = seqan::row(align,1);
	return length(source(row1))- countTotalGaps(row2);
}

template <typename TAlign>
unsigned getInsertSize(TAlign& align) noexcept
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	const TRow &row1 = seqan::row(align,0);
	const TRow &row2 = seqan::row(align,1);
	const unsigned seq1_length = length(source(row1));
	const unsigned seq2_length = length(source(row2));
	// Calculate overlap and overhangs.
	const unsigned overlap = seq1_length - countTotalGaps(row2);
	const unsigned seq2l = seqan::countGaps(seqan::begin(row1)); // Overhang of sequence 2 = Gaps at start of row 1.
	const unsigned seq1r = countTotalGaps(row2) - seqan::countGaps(seqan::begin(row2)); // Overhang of sequence 1 = Gaps at end of row 2.
	// Insert size: Add sequence lengths, subtract common region (overlap)
	// and subtract overhangs left and right of insert.
	return seq1_length + seq2_length - overlap - (seq2l + seq1r);
}

namespace AlignAlgorithm
{
    struct NeedlemanWunsch {};
    struct Menkuec {};
}

template <typename TSeq, typename TAdapter>
void alignPair(std::pair<int, seqan::Align<TSeq> >& ret, const TSeq& seq1, const TAdapter& seq2, 
    const int leftOverhang, const int rightOverhang, const AlignAlgorithm::NeedlemanWunsch&) noexcept
{
    seqan::resize(rows(ret.second), 2);
    seqan::assignSource(row(ret.second, 0), seq1);
    seqan::assignSource(row(ret.second, 1), seq2);

    // dont allow gaps by setting the gap penalty high
    const TScore adapterScore(-100);

    // banded alignment
    const int shiftStartPos = -leftOverhang;
    const int shiftEndPos = length(seq1) - length(seq2) + rightOverhang;
    if (shiftEndPos < shiftStartPos)
    {
        ret.first = -1;
        return;
    }
    // allow gaps on all corners, because we use banded alignment to constraint the search space
    seqan::AlignConfig<true, true, true, true> config; 
    ret.first = globalAlignment(ret.second, adapterScore, config, shiftStartPos, shiftEndPos, seqan::LinearGaps());
}

/*
- shifts adapterTemplate against sequence
- calculate score for each shift position
  - +1 for same base, -1 for mismatch, +0 for N
- return score of the shift position, where the errorRate was minimal
*/
template <typename TSeq, typename TAdapter>
void alignPair(std::pair<int, seqan::Align<TSeq> >& ret, const TSeq& seq1, const TAdapter& seq2, 
        const int leftOverhang, const int rightOverhang, const AlignAlgorithm::Menkuec&) noexcept
{
    seqan::resize(rows(ret.second), 2);
    seqan::assignSource(row(ret.second, 0), seq1);
    seqan::assignSource(row(ret.second, 1), seq2);

    const int shiftStartPos = -leftOverhang;
    const int shiftEndPos = length(seq1) - length(seq2) + rightOverhang;
    const auto lenSeq1 = length(seq1);
    const auto lenSeq2 = length(seq2);
    int shiftPos = shiftStartPos;
    int bestShiftPos = shiftStartPos;
    int bestScore = std::numeric_limits<int>::min();
    float bestErrorRate = std::numeric_limits<float>::max();
    unsigned int bestOverlap = 0;

    if (shiftEndPos < shiftStartPos)
    {
        ret.first = bestScore; // invalid constraints
        return;
    }
    while (shiftPos <= shiftEndPos)
    {
        const unsigned int overlapNegativeShift = std::min(shiftPos + lenSeq2, lenSeq1);
        const unsigned int overlapPositiveShift = std::min(lenSeq1 - shiftPos, lenSeq2);
        const unsigned int overlap = std::min(overlapNegativeShift, overlapPositiveShift);
        const unsigned int overlapStart = std::max(shiftPos, 0);
        int score = 0;
        for (unsigned int pos = 0; pos < overlap; ++pos)
        {
            if (seq2[pos + std::min(0,shiftPos)*(-1)] == seq1[overlapStart + pos])
                ++score;
            else if (seq1[overlapStart + pos] != 'N')
                --score;
        }
        const float errorRate = static_cast<float>((overlap-score)/2) / static_cast<float>(overlap);
        if (errorRate < bestErrorRate || (errorRate == bestErrorRate && overlap > bestOverlap))
        {
            bestShiftPos = shiftPos;
            bestScore = score;
            bestErrorRate = errorRate;
            bestOverlap = overlap;
        }
        ++shiftPos;
    }
    if (bestShiftPos < 0)
    {
        seqan::insertGaps(row(ret.second, 0), 0, - bestShiftPos); // top left
        seqan::insertGaps(row(ret.second, 0), lenSeq1 - bestShiftPos, std::max<int>(0, lenSeq2 - lenSeq1 + bestShiftPos)); // top right
        seqan::insertGaps(row(ret.second, 1), lenSeq2, std::max<int>(0,lenSeq1 - lenSeq2 - bestShiftPos)); // bottom right
    }
    else
    {
        seqan::insertGaps(row(ret.second, 0), lenSeq1, std::max<int>(0, lenSeq2 - lenSeq1 + bestShiftPos)); // top right
        seqan::insertGaps(row(ret.second, 1), 0, bestShiftPos); // bottom left
        seqan::insertGaps(row(ret.second, 1), lenSeq2+bestShiftPos, std::max<int>(0,lenSeq1 - lenSeq2 - bestShiftPos)); // bottom right
    }
    ret.first = bestScore;
}

// used only for testing and paired end data
template <typename TSeq, typename TAdapter>
void alignPair(std::pair<int, seqan::Align<TSeq> >& ret, const TSeq& seq1, const TAdapter& seq2, const AlignAlgorithm::NeedlemanWunsch&) noexcept
{
    seqan::resize(rows(ret.second), 2);
    seqan::assignSource(row(ret.second, 0), seq1);
    seqan::assignSource(row(ret.second, 1), seq2);

    // dont allow gaps by setting the gap penalty high
    const TScore adapterScore(-100);

	// global alignment, not banded
    seqan::AlignConfig<true, true, true, true> config;
    ret.first = globalAlignment(ret.second, adapterScore, config, seqan::NeedlemanWunsch());
}

template <typename TSeq>
unsigned stripPair(TSeq& seq1, TSeq& seq2) noexcept
{
    // When aligning the two sequences, the complementary sequence is reversed and
    // complemented, so we have an overlap alignment with complementary bases being the same.
    using TAlphabet = typename seqan::Value<TSeq>::Type;
    using TReverseComplement = typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type;
    TReverseComplement mod(seq2);
    typedef seqan::Align<TSeq> TAlign;
    std::pair<int, TAlign> ret;
    alignPair(ret, seq1, mod, AlignAlgorithm::NeedlemanWunsch());
    const auto& score = ret.first;
    const TAlign& align = ret.second;
    // Use the overlap of the two sequences to determine the end position.
    const unsigned overlap = getOverlap(align);
    const unsigned mismatches = (overlap - score) / 2;
    // We require a certain correct overlap to exclude spurious hits.
    // (Especially reverse 3'->5' alignments not caused by adapters.)
    if (overlap <= 5 || mismatches > overlap * 0.15)
    {
        return 0;
    }
    // Get actual size of the insert (possible to determine from overlap etc.).
    const unsigned insert = getInsertSize(align);
    // Now cut both sequences to insert size (no cuts happen if they are smaller)
    if (length(seq1) > insert)
    {
        seqan::erase(seq1, insert, length(seq1));
    }
    if (length(seq2) > insert)
    {
        seqan::erase(seq2, insert, length(seq2));
    }
    return insert;
}

//Version for automatic matching options
inline bool isMatch(const unsigned int overlap, const int mismatches, const AdapterMatchSettings &adatperMatchSettings) noexcept
{
    if (overlap == 0)
        return false;
    if (adatperMatchSettings.errorRate > 0)
        return overlap >= adatperMatchSettings.min_length && (static_cast<double>(mismatches) / static_cast<double>(overlap)) <= adatperMatchSettings.errorRate;
    else
        return overlap >= adatperMatchSettings.min_length && mismatches <= adatperMatchSettings.errors;
}

enum adapterDirection : bool
{
    reverse, 
    forward
};

template<bool _val>
struct TagAdapter
{
    static const bool value = _val;
};

template<bool _direction>
struct StripAdapterDirection
{
    static const bool value = _direction;
};

template <typename TSeq, typename TAdapters, typename TStripAdapterDirection>
unsigned stripAdapter(TSeq& seq, AdapterTrimmingStats& stats, TAdapters const& adapters, AdapterMatchSettings const& spec,
    const TStripAdapterDirection&)
{
    using TAlign = seqan::Align<TSeq>;
    AlignAlgorithm::Menkuec alignAlgorithm;

    unsigned removed{ 0 };
    std::vector<std::tuple<unsigned int, unsigned int, TAlign, AdapterItem>> matches;
    matches.reserve(spec.times);
    for (unsigned int n = 0;n < spec.times; ++n)
    {
        matches.clear();
        {
            std::pair<int, TAlign> ret;
            for (auto const& adapterItem : adapters)
            {
                //std::cout << "seq: " << seq << std::endl;
                //std::cout << "adapter: " << adapterItem.seq << std::endl;
                if (static_cast<unsigned>(length(adapterItem.seq)) < spec.min_length)
                    continue;
                if (TStripAdapterDirection::value == adapterDirection::reverse && adapterItem.reverse == false)
                    continue;
                if (TStripAdapterDirection::value == adapterDirection::forward && adapterItem.reverse == true)
                    continue;

                const auto& adapterSequence = adapterItem.seq;
                const unsigned int oppositeEndOverhang = adapterItem.anchored == true ? length(adapterSequence) - length(seq) : adapterItem.overhang;
                const unsigned int sameEndOverhang = adapterItem.anchored == true ? 0 : length(adapterItem.seq) - spec.min_length;
                if (adapterItem.adapterEnd == AdapterItem::end3)
                    alignPair(ret, seq, adapterSequence, oppositeEndOverhang, sameEndOverhang, alignAlgorithm);
                else
                    alignPair(ret, seq, adapterSequence, sameEndOverhang, oppositeEndOverhang, alignAlgorithm);

                const int score = ret.first;
                if (score < 0)
                    continue;
                const unsigned int overlap = getOverlap(ret.second);
                const int mismatches = (overlap - score) / 2;
                //std::cout << "score: " << ret.first << " er: " << (float)mismatches/overlap << " overlap: " << overlap << " mismatches: " << mismatches << std::endl;
                //std::cout << ret.second << std::endl;

                if (isMatch(overlap, mismatches, spec))
                    matches.push_back(std::make_tuple(score, overlap, ret.second, std::ref(adapterItem)));
            }
        }
        if (matches.empty())
            return removed;

        // select best matching adapter
        auto maxIt = matches.cbegin();
        for (auto it = matches.cbegin(); it != matches.cend();++it)
            if (std::get<0>(*it)>std::get<0>(*maxIt))
                maxIt = it;

        // erase best matching adapter from sequence
        const auto score = std::get<0>(*maxIt);
        const auto overlap = std::get<1>(*maxIt);
        const TAlign align = std::get<2>(*maxIt);
        AdapterItem adapterItem = std::get<3>(*maxIt);
        const auto mismatches = (overlap - score) / 2;
        unsigned eraseStart = 0;
        unsigned eraseEnd = 0;
        if (adapterItem.adapterEnd == AdapterItem::end3)
        {
            eraseStart = toViewPosition(row(align, 1), 0);
            eraseEnd = length(seq);
            //std::cout << "adapter start position: " << toViewPosition(row(ret.second, 1), 0) << std::endl;
        }
        else
        {
            eraseStart = 0;
            eraseEnd = std::min<unsigned>(length(seq),toViewPosition(row(align, 1), length(adapterItem.seq)) - toViewPosition(row(align, 0), 0));
            //std::cout << "adapter end position: " << toViewPosition(row(ret.second, 1), length(adapterItem.seq)) - toViewPosition(row(ret.second, 0), 0) << std::endl;
        }

        //std::cout << "unstripped seq: " << seq << " erase: " << eraseStart << " " << eraseEnd;
        seqan::erase(seq, eraseStart, eraseEnd);
        removed += eraseEnd - eraseStart;
        //std::cout << "\nstripped seq  : " << seq << std::endl;

        // update statistics

        const auto statisticLen = eraseEnd - eraseStart;
        if (stats.removedLength.size() < statisticLen)
            stats.removedLength.resize(statisticLen);
        if (stats.removedLength[statisticLen - 1].size() < mismatches + 1)
            stats.removedLength[statisticLen - 1].resize(mismatches + 1);
        ++stats.removedLength[statisticLen - 1][mismatches];

        if (stats.numRemoved.size() < adapterItem.id + 1)
        {
            std::cout << "error: numRemoved too small!" << std::endl;
            throw(std::runtime_error("error: numRemoved too small!"));
        }
        ++stats.numRemoved[adapterItem.id];

        stats.overlapSum += overlap;
        stats.maxOverlap = std::max(stats.maxOverlap, overlap);
        stats.minOverlap = std::min(stats.minOverlap, overlap);
        // dont try more adapter trimming if the read is too short already
        if (static_cast<unsigned>(length(seq)) < spec.min_length)
            return removed;
    }
    return removed;
}

template < template <typename> class TRead, typename TSeq, typename TAdaptersArray, typename TSpec, typename TTagAdapter,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value> >
void stripAdapterBatch(std::vector<TRead<TSeq>>& reads, TAdaptersArray const& adapters, TSpec const& spec, const bool pairedNoAdapterFile,
    AdapterTrimmingStats& stats, TTagAdapter, bool = false) noexcept(!TTagAdapter::value)
{
    (void)pairedNoAdapterFile;
    for(auto& read: reads)
    {
        if (seqan::empty(read.seq))
            continue;
        const unsigned over = stripAdapter(read.seq, stats, adapters, spec, StripAdapterDirection<adapterDirection::forward>());
        if (TTagAdapter::value && over != 0)
            insertAfterFirstToken(read.id, ":AdapterRemoved");
    }
    return;
}

// pairedEnd adapters will be trimmed in single mode, each seperately
template < template <typename> class TRead, typename TSeq, typename TAdaptersArray, typename TSpec, typename TTagAdapter,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value> >
    void stripAdapterBatch(std::vector<TRead<TSeq>>& reads, TAdaptersArray const& adapters, TSpec const& spec, const bool pairedNoAdapterFile,
        AdapterTrimmingStats& stats, TTagAdapter) noexcept(!TTagAdapter::value)
{
    for (auto& read : reads)
    {
        if (seqan::empty(read.seq))
            continue;
        unsigned over = 0;
        if (pairedNoAdapterFile)
        {
            stripPair(read.seq, read.seqRev);
        }
        else
        {
            over = stripAdapter(read.seq, stats, adapters, spec, StripAdapterDirection<adapterDirection::forward>());
            if (!seqan::empty(read.seqRev))
                over += stripAdapter(read.seqRev, stats, adapters, spec, StripAdapterDirection<adapterDirection::reverse>());
        }
        if (TTagAdapter::value && over != 0)
            insertAfterFirstToken(read.id, ":AdapterRemoved");
    }
    return;
}

