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
#include <xmmintrin.h>
#include <nmmintrin.h>

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
namespace seqan {
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

using TAdapterSequence = std::string;

struct AdapterItem;
typedef std::vector< AdapterItem > AdapterSet;

std::string ReverseComplement(std::string str)
{
    // TODO implement function
    return str;
}

struct AdapterItem
{
    enum AdapterEnd
    {
        end3,
        end5,
    };
    AdapterEnd adapterEnd;
    unsigned char overhang;
    unsigned char id;
    bool anchored;
    bool reverse;

    AdapterItem() : adapterEnd(end3), overhang(0), id(0), anchored(false), reverse(false), len(0) {};
    AdapterItem(const TAdapterSequence &adapter) : adapterEnd(end3), overhang(0), id(0), anchored(false), reverse(false), seq(adapter), len(length(adapter)) {};
    AdapterItem(const TAdapterSequence &adapter, const AdapterEnd adapterEnd, const unsigned overhang, const unsigned id, const bool anchored, const bool reverse)
        : adapterEnd(adapterEnd), overhang(overhang), id(id), anchored(anchored), reverse(reverse), seq(adapter), len(length(adapter)) {};

    void setSeq(TAdapterSequence newSeq) noexcept
    { 
        seq = newSeq; 
        len = length(seq); 
    };
    const TAdapterSequence& getSeq() const noexcept
    {
        return seq;
    };
    unsigned char getLen() const noexcept
    {
        return len;
    };
    
    AdapterItem getReverseComplement() const noexcept
    {
        auto seqCopy = seq;
        return AdapterItem(ReverseComplement(seqCopy), adapterEnd, overhang, id, anchored, reverse);
    }
private:
    TAdapterSequence seq;
    unsigned char len;
};

// Define scoring function type.
typedef seqan::Score<int, seqan::ScoreMatrix<seqan::Dna5, seqan::AdapterScoringMatrix> > TScore;

struct AdapterMatchSettings
{
    AdapterMatchSettings(const char m, const int e, const double er, const unsigned char oh, const unsigned char times) : min_length(m), errors(e), overhang(oh), errorRate(er), times(times)
    {}
    AdapterMatchSettings() : min_length(0), errors(0), overhang(0), errorRate(0), times(1) {};

    unsigned char min_length; //The minimum length of the overlap.
    unsigned int errors;     //The maximum number of errors we allow.
    unsigned char overhang;
    double errorRate;  //The maximum number of errors allowed per overlap
    unsigned char times;
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
    const TRow &row1 = seqan::row(align, 0);
    const TRow &row2 = seqan::row(align, 1);
    return length(source(row1)) - countTotalGaps(row2);
}

template <typename _TReadLen>
struct AlignResult
{
    static const typename std::make_signed<_TReadLen>::type noMatch = std::numeric_limits<typename std::make_signed<_TReadLen>::type>::min();

    AlignResult() : shiftPos(0), score(std::numeric_limits<typename std::make_signed<_TReadLen>::type>::min()), matches(0), ambiguous(0), overlap(0), errorRate(1) {};
    typename std::make_signed<_TReadLen>::type shiftPos;
    typename std::make_signed<_TReadLen>::type score;
    typename std::make_unsigned<_TReadLen>::type matches;
    typename std::make_unsigned<_TReadLen>::type mismatches;
    typename std::make_unsigned<_TReadLen>::type ambiguous;
    typename std::make_unsigned<_TReadLen>::type overlap;
    float errorRate;
};

template <typename TAlign>
unsigned getInsertSize(TAlign& align) noexcept
{
    typedef typename seqan::Row<TAlign>::Type TRow;
    const TRow &row1 = seqan::row(align, 0);
    const TRow &row2 = seqan::row(align, 1);
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


template <typename TSeq, typename TAdapter, typename TAlignResult>
void alignPair(TAlignResult &res, const TSeq& seq1, const TAdapter& seq2,
    const int leftOverhang, const int rightOverhang, const AlignAlgorithm::NeedlemanWunsch&) noexcept
{
    seqan::Align<TSeq> align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq1);
    seqan::assignSource(row(align, 1), seq2);

    // dont allow gaps by setting the gap penalty high
    const TScore adapterScore(-100);

    // banded alignment
    const int shiftStartPos = -leftOverhang;
    const int shiftEndPos = length(seq1) - length(seq2) + rightOverhang;
    if (shiftEndPos < shiftStartPos)
    {
        res.score = TAlignResult::noMatch;
        return;
    }
    // allow gaps on all corners, because we use banded alignment to constraint the search space
    seqan::AlignConfig<true, true, true, true> config;
    res.score = globalAlignment(align, adapterScore, config, shiftStartPos, shiftEndPos, seqan::LinearGaps());
    res.shiftPos = toViewPosition(row(align, 1), 0);
    res.overlap = getOverlap(align);
    res.ambiguous = 0;  // TODO: calculate ambiguous matchings
    res.matches = (res.score + res.overlap)/2;
    res.errorRate = static_cast<float>(res.overlap - res.matches) / static_cast<float>(res.overlap);
}

const __m128i ONE_8 = _mm_set1_epi64x(0x01);
const __m128i ONE_16 = _mm_set1_epi64x(0x0101);
const __m128i ONE_24 = _mm_set1_epi64x(0x010101);
const __m128i ONE_32 = _mm_set1_epi64x(0x01010101);
const __m128i ONE_40 = _mm_set1_epi64x(0x0101010101);
const __m128i ONE_48 = _mm_set1_epi64x(0x010101010101);
const __m128i ONE_56 = _mm_set1_epi64x(0x01010101010101);
const __m128i ONE_128 = _mm_set1_epi8(1);
const __m256i ONE_256 = _mm256_set1_epi8(1);

const __m128i ZERO_128 = _mm_set1_epi8(0);
const __m256i ZERO_256 = _mm256_set1_epi8(0);
//const __m128i N_128 = _mm_set1_epi8(0x04);
const __m128i N_128 = _mm_set1_epi8('N');
const __m256i N_256 = _mm256_set1_epi8('N');

// vector access to SSE registers is a microsoft specialty
#ifdef _MSC_VER
    #define VECTOR_ACCESS
#endif

inline size_t popcnt64(__m128i value) noexcept
{
    //value = _mm_sad_epu8(ZERO_128, value);
    //return _mm_extract_epi16(value, 0);
#ifdef VECTOR_ACCESS
    return _mm_popcnt_u64(value.m128i_u64[0]);
#else
    return _mm_popcnt_u64(_mm_extract_epi64(value,0));
#endif
}

inline size_t popcnt128(__m128i value) noexcept
{
//    value = _mm_sad_epu8(ZERO_128, value);
//    return _mm_extract_epi16(value, 0) + _mm_extract_epi16(value, 4);
#ifdef VECTOR_ACCESS
    return _mm_popcnt_u64(value.m128i_u64[0]) + _mm_popcnt_u64(value.m128i_u64[1]);
#else
    return _mm_popcnt_u64(_mm_extract_epi64(value, 0)) + _mm_popcnt_u64(_mm_extract_epi64(value, 1));
#endif
}

inline size_t popcnt256(__m256i value) noexcept
{
#ifdef VECTOR_ACCESS
    return _mm_popcnt_u64(value.m256i_u64[0]) +
        _mm_popcnt_u64(value.m256i_u64[1]) +
        _mm_popcnt_u64(value.m256i_u64[2]) +
        _mm_popcnt_u64(value.m256i_u64[3]);
#else
    return _mm_popcnt_u64(_mm256_extract_epi64(value, 0)) +
        _mm_popcnt_u64(_mm256_extract_epi64(value, 1)) +
        _mm_popcnt_u64(_mm256_extract_epi64(value, 2)) +
        _mm_popcnt_u64(_mm256_extract_epi64(value, 3));
#endif
}

template <unsigned int N>
struct compareAdapter
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        // SSE2 code
        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_8));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_8));
        
        // SSE4.2 code 
        //ambiguous += _mm_popcnt_u32(_mm_extract_epi8(_mm_and_si128(NMask, ONE_128), 0));
        //matches += _mm_popcnt_u32(_mm_extract_epi8(_mm_and_si128(matchesMask, ONE_128), 0));

        ++adapterIterator;
        ++readIterator;
        compareAdapter<N-1>::apply(readIterator, adapterIterator, matches, ambiguous);
    }
};

template <>
struct compareAdapter<0>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        (void)readIterator;
        (void)adapterIterator;
        (void)matches;
        (void)ambiguous;
    }
};

template <>
struct compareAdapter<2>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_16));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_16));
        //ambiguous += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(NMask, ONE_128),0));
        //matches += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(matchesMask, ONE_128), 0));
        readIterator += 2;
        adapterIterator += 2;
    }
};

template <>
struct compareAdapter<3>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_24));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_24));
        //ambiguous += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(NMask, ONE_128),0));
        //matches += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(matchesMask, ONE_128), 0));
        readIterator += 3;
        adapterIterator += 3;
    }
};

template <>
struct compareAdapter<4>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_32));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_32));
        //ambiguous += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(NMask, ONE_128),0));
        //matches += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(matchesMask, ONE_128), 0));
        readIterator += 4;
        adapterIterator += 4;
    }
};

template <>
struct compareAdapter<5>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_40));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_40));
        //ambiguous += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(NMask, ONE_128),0));
        //matches += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(matchesMask, ONE_128), 0));
        readIterator += 5;
        adapterIterator += 5;
    }
};

template <>
struct compareAdapter<6>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_48));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_48));
        //ambiguous += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(NMask, ONE_128),0));
        //matches += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(matchesMask, ONE_128), 0));
        readIterator += 6;
        adapterIterator += 6;
    }
};

template <>
struct compareAdapter<7>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_56));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_56));
        //ambiguous += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(NMask, ONE_128),0));
        //matches += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(matchesMask, ONE_128), 0));
        readIterator += 7;
        adapterIterator += 7;
    }
};

template <>
struct compareAdapter<8>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128), _mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt64(_mm_and_si128(NMask, ONE_128));
        matches += popcnt64(_mm_and_si128(matchesMask, ONE_128));
        //ambiguous += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(NMask, ONE_128),0));
        //matches += _mm_popcnt_u64(_mm_extract_epi64(_mm_and_si128(matchesMask, ONE_128), 0));
        readIterator += 8;
        adapterIterator += 8;
    }
};

template <>
struct compareAdapter<16>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m128i read = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*readIterator)));
        const __m128i adapter = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&(*adapterIterator)));

        const __m128i NMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, N_128),_mm_cmpeq_epi8(adapter, N_128)));
        const __m128i matchesMask = _mm_sub_epi8(ZERO_128, _mm_or_si128(_mm_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt128(_mm_and_si128(NMask, ONE_128));
        matches += popcnt128(_mm_and_si128(matchesMask, ONE_128));
        readIterator += 16;
        adapterIterator += 16;
    }
};

template <>
struct compareAdapter<32>
{
    template <typename TReadIterator, typename TAdapterIterator, typename TCounter>
    inline static void apply(TReadIterator& readIterator, TAdapterIterator& adapterIterator, TCounter& matches, TCounter& ambiguous) noexcept
    {
        const __m256i read = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&(*readIterator)));
        const __m256i adapter = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&(*adapterIterator)));

        const __m256i NMask = _mm256_sub_epi8(ZERO_256, _mm256_or_si256(_mm256_cmpeq_epi8(read, N_256), _mm256_cmpeq_epi8(adapter, N_256)));
        const __m256i matchesMask = _mm256_sub_epi8(ZERO_256, _mm256_or_si256(_mm256_cmpeq_epi8(read, adapter), NMask));

        ambiguous += popcnt256(_mm256_and_si256(NMask, ONE_256));
        matches += popcnt256(_mm256_and_si256(matchesMask, ONE_256));
        readIterator += 32;
        adapterIterator += 32;
    }
};

/*
- shifts adapterTemplate against sequence
- calculate score for each shift position
- +1 for same base, -1 for mismatch, +0 for N
- return score of the shift position, where the errorRate was minimal
*/
template <typename TSeq, typename TAdapter, typename TAlignResult>
void alignPair(TAlignResult& ret, const TSeq& read, const TAdapter& adapter,
    const int leftOverhang, const int rightOverhang, const AlignAlgorithm::Menkuec&) noexcept
{
    const auto lenRead = length(read);
    const auto lenAdapter = length(adapter);
    const int shiftStartPos = -leftOverhang;
    const int shiftEndPos = lenRead - lenAdapter + rightOverhang;
    int shiftPos = shiftStartPos;

    ret = TAlignResult();
    const std::string::const_iterator readBeginIterator = read.begin();
    const std::string::const_iterator adapterBeginIterator = adapter.begin();
    std::string::const_iterator readIterator = readBeginIterator;
    std::string::const_iterator adapterIterator = adapterBeginIterator;
    while (shiftPos <= shiftEndPos)
    {
        const unsigned int overlapNegativeShift = std::min<unsigned int>(shiftPos + lenAdapter, lenRead);
        const unsigned int overlapPositiveShift = std::min<unsigned int>(lenRead - shiftPos, lenAdapter);
        const unsigned int overlap = std::min<unsigned int>(overlapNegativeShift, overlapPositiveShift);
        const unsigned int overlapStart = std::max<int>(shiftPos, 0);
        unsigned int matches = 0;
        unsigned int ambiguous = 0;
        unsigned int remaining = overlap;
        readIterator = readBeginIterator + overlapStart;
        adapterIterator = adapterBeginIterator + std::min(0, shiftPos)*(-1);

        while (remaining >= 32)
        {
            compareAdapter<32>::apply(readIterator, adapterIterator, matches, ambiguous);
            remaining -= 32;
        }
        if (remaining >= 16)
        {
            compareAdapter<16>::apply(readIterator, adapterIterator, matches, ambiguous);
            remaining -= 16;
        }
        if (remaining >= 8)
        {
            compareAdapter<8>::apply(readIterator, adapterIterator, matches, ambiguous);
            remaining -= 8;
        }
        switch (remaining)
        {
        case 0:
            break;
        case 1:
            compareAdapter<1>::apply(readIterator, adapterIterator, matches, ambiguous);
            break;
        case 2:
            compareAdapter<2>::apply(readIterator, adapterIterator, matches, ambiguous);
            break;
        case 3:
            compareAdapter<3>::apply(readIterator, adapterIterator, matches, ambiguous);
            break;
        case 4:
            compareAdapter<4>::apply(readIterator, adapterIterator, matches, ambiguous);
            break;
        case 5:
            compareAdapter<5>::apply(readIterator, adapterIterator, matches, ambiguous);
            break;
        case 6:
            compareAdapter<6>::apply(readIterator, adapterIterator, matches, ambiguous);
            break;
        case 7:
            compareAdapter<7>::apply(readIterator, adapterIterator, matches, ambiguous);
            break;
        }
        const float errorRate = static_cast<float>(overlap - matches - ambiguous) / static_cast<float>(overlap);
        if (errorRate < ret.errorRate || (errorRate == ret.errorRate && overlap > ret.overlap))
        {
            ret.matches = matches;
            ret.ambiguous = ambiguous;
            ret.errorRate = errorRate;
            ret.shiftPos = shiftPos;
            ret.overlap = overlap;
        }
        ++shiftPos;
    }
    if (ret.matches != 0)
    {
        ret.mismatches = ret.overlap - ret.matches - ret.ambiguous;
        ret.score = 2 * ret.matches - ret.overlap + ret.ambiguous;
    }
}



template <typename TSeq, typename TAdapter, typename TAlignResult>
void alignPair(TAlignResult &res, const TSeq& seq1, const TAdapter& seq2, const AlignAlgorithm::NeedlemanWunsch&) noexcept
{
    seqan::Align<TSeq> align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq1);
    seqan::assignSource(row(align, 1), seq2);

    // dont allow gaps by setting the gap penalty high
    const TScore adapterScore(-100);

    // allow gaps on all corners, because we use banded alignment to constraint the search space
    seqan::AlignConfig<true, true, true, true> config;
    res.score = globalAlignment(align, adapterScore, config, seqan::LinearGaps());
    res.shiftPos = toViewPosition(row(align, 1), 0);
    res.overlap = getOverlap(align);
    res.ambiguous = 0;  // TODO: calculate ambiguous matchings
    res.matches = (res.score + res.overlap) / 2;
    res.errorRate = static_cast<float>(res.overlap - res.matches) / static_cast<float>(res.overlap);
}

template <typename TSeq>
unsigned stripPair(TSeq& seq1, TSeq& seq2) noexcept
{
    // When aligning the two sequences, the complementary sequence is reversed and
    // complemented, so we have an overlap alignment with complementary bases being the same.
    using TAlphabet = typename seqan::Value<TSeq>::Type;
    using TReverseComplement = typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type;
    TReverseComplement mod(seq2);
    AlignResult<unsigned int> ret;
    alignPair(ret, seq1, mod, AlignAlgorithm::NeedlemanWunsch());
    const auto& score = ret.score;
    // Use the overlap of the two sequences to determine the end position.
    const unsigned overlap = ret.overlap;
    const unsigned mismatches = (overlap - score) / 2;
    // We require a certain correct overlap to exclude spurious hits.
    // (Especially reverse 3'->5' alignments not caused by adapters.)
    if (overlap <= 5 || mismatches > overlap * 0.15)
    {
        return 0;
    }
    // Get actual size of the insert (possible to determine from overlap etc.).
    const unsigned insert = length(seq1) + length(seq1) - overlap;   // TODO: check gaps in getInsertSize(align);
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
inline bool isMatch(const unsigned int overlap, const unsigned int mismatches, const AdapterMatchSettings &adatperMatchSettings) noexcept
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

void Dna5ToStdString(std::string& dest, const seqan::Dna5QString &source) noexcept
{
    auto len = length(source);
    dest.resize(len);
    for(size_t n = 0; n < len; n++)
    {
        switch (((unsigned char)source[n]) & 0x07)
        {
        case 0:
            dest[n] = 'A';
            break;
        case 1:
            dest[n] = 'C';
            break;
        case 2:
            dest[n] = 'G';
            break;
        case 3:
            dest[n] = 'T';
            break;
        case 4:
            dest[n] = 'N';
            break;
        }
    }
}

template <typename TStats>
struct TlsBlockAdapterTrimming
{
    TlsBlockAdapterTrimming(TStats& stats, const AdapterTrimmingParams& params) : stats(stats), params(params) {};

    TStats& stats;
    const AdapterTrimmingParams& params; // can use ref here, bcs read only does not cause false sharing
    std::string tlsString;
};

// convenience wrapper
template <typename TSeq, typename TAdapters, typename TReadLen, typename TStripAdapterDirection>
unsigned stripAdapter(TSeq& seq, AdapterTrimmingStats<TReadLen>& stats, TAdapters const& adapters, AdapterMatchSettings const& spec,
    const TStripAdapterDirection& stripDirection)
{
    AdapterTrimmingParams params;
    params.adapters = adapters;
    params.mode = spec;
    TlsBlockAdapterTrimming<AdapterTrimmingStats<TReadLen>> tlsBlock(stats, params);
    return stripAdapter(seq, tlsBlock, stripDirection);
}

template <typename TSeq, typename TStripAdapterDirection, typename TlsBlock>
unsigned stripAdapter(TSeq& seq, TlsBlock& tlsBlock, const TStripAdapterDirection&)
{
    AlignAlgorithm::Menkuec alignAlgorithm;

    using TReadLen = decltype(tlsBlock.stats.overlapSum);
    unsigned removedTotal{ 0 };
    AlignResult<TReadLen> alignResult;  // small object, created on stack
    unsigned removedTotalOld = 0;
    TReadLen lenSeq = length(seq);

    Dna5ToStdString(tlsBlock.tlsString, seq);

    for (unsigned int n = 0;n < tlsBlock.params.mode.times; ++n)
    {
        alignResult.score = AlignResult<TReadLen>::noMatch;
        for (auto const& adapterItem : tlsBlock.params.adapters)
        {
            //if (static_cast<unsigned>(length(adapterItem.seq)) < spec.min_length)
              //  continue;
            if ((TStripAdapterDirection::value == adapterDirection::reverse && adapterItem.reverse == false) ||
                (TStripAdapterDirection::value == adapterDirection::forward && adapterItem.reverse == true))
                continue;

            const auto& adapterSequence = adapterItem.getSeq();
            const auto lenAdapter = adapterItem.getLen();

            const int oppositeEndOverhang = adapterItem.anchored == true ? lenAdapter - lenSeq : adapterItem.overhang;
            const int sameEndOverhang = adapterItem.anchored == true ? 0 : lenAdapter - tlsBlock.params.mode.min_length;
            if (adapterItem.adapterEnd == AdapterItem::end3)
                alignPair(alignResult, tlsBlock.tlsString, adapterSequence, oppositeEndOverhang, sameEndOverhang, alignAlgorithm);
            else
                alignPair(alignResult, tlsBlock.tlsString, adapterSequence, sameEndOverhang, oppositeEndOverhang, alignAlgorithm);

            if (isMatch(alignResult.overlap, alignResult.mismatches, tlsBlock.params.mode))
            {
                TReadLen eraseStart = 0;
                TReadLen eraseEnd = 0;
                if (adapterItem.adapterEnd == AdapterItem::end3)
                {
                    eraseStart = alignResult.shiftPos;
                    eraseEnd = lenSeq;
                }
                else
                {
                    eraseStart = 0;
                    eraseEnd = std::min<TReadLen>(lenSeq, alignResult.shiftPos + lenAdapter);
                }

                seqan::erase(seq, eraseStart, eraseEnd);
                tlsBlock.tlsString.erase(eraseStart, eraseEnd);
                TReadLen removed = eraseEnd - eraseStart;
                removedTotal += removed;
                lenSeq -= removed;

                // update statistics        
                const auto statisticLen = removed;
                if (tlsBlock.stats.removedLength.size() < statisticLen)
                    tlsBlock.stats.removedLength.resize(statisticLen);
                if (tlsBlock.stats.removedLength[statisticLen - 1].size() < static_cast<size_t>(alignResult.mismatches + 1))
                    tlsBlock.stats.removedLength[statisticLen - 1].resize(alignResult.mismatches + 1);
                ++tlsBlock.stats.removedLength[statisticLen - 1][alignResult.mismatches];

                if (tlsBlock.stats.numRemoved.size() < static_cast<size_t>(adapterItem.id + 1))
                {
                    std::cout << "error: numRemoved too small!" << std::endl;
                    throw(std::runtime_error("error: numRemoved too small!"));
                }
                ++tlsBlock.stats.numRemoved[adapterItem.id];

                tlsBlock.stats.overlapSum += alignResult.overlap;
                tlsBlock.stats.maxOverlap = std::max(tlsBlock.stats.maxOverlap, alignResult.overlap);
                tlsBlock.stats.minOverlap = std::min(tlsBlock.stats.minOverlap, alignResult.overlap);
            }
        }
        if (removedTotal == removedTotalOld)
            return removedTotal;
        removedTotalOld = removedTotal;

        // dont try more adapter trimming if the read is too short already
        if (static_cast<TReadLen>(lenSeq) < tlsBlock.params.mode.min_length)
            return removedTotal;
    }
    return removedTotal;
}

template < template <typename> class TRead, typename TSeq, typename TlsBlock, typename TTagAdapter,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value> >
    void stripAdapterBatch(std::vector<TRead<TSeq>>& reads, TlsBlock& tlsBlock, TTagAdapter, bool = false) noexcept(!TTagAdapter::value)
{
    for (auto& read : reads)
    {
        if (seqan::empty(read.seq))
            continue;
        const unsigned over = stripAdapter(read.seq, tlsBlock, StripAdapterDirection<adapterDirection::forward>());
        if (TTagAdapter::value && over != 0)
            insertAfterFirstToken(read.id, ":AdapterRemoved");
    }
    return;
}

// pairedEnd adapters will be trimmed in single mode, each seperately
template < template <typename> class TRead, typename TSeq, typename TlsBlock, typename TTagAdapter,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value> >
    void stripAdapterBatch(std::vector<TRead<TSeq>>& reads, TlsBlock& tlsBlock, TTagAdapter) noexcept(!TTagAdapter::value)
{
    for (auto& read : reads)
    {
        if (seqan::empty(read.seq))
            continue;
        unsigned over = 0;
        if (tlsBlock.params.pairedNoAdapterFile)
        {
            stripPair(read.seq, read.seqRev);
        }
        else
        {
            over = stripAdapter(read.seq, tlsBlock, StripAdapterDirection<adapterDirection::forward>());
            if (!seqan::empty(read.seqRev))
                over += stripAdapter(read.seqRev, tlsBlock, StripAdapterDirection<adapterDirection::reverse>());
        }
        if (TTagAdapter::value && over != 0)
            insertAfterFirstToken(read.id, ":AdapterRemoved");
    }
    return;
}

