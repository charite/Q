#ifndef BAM_RECORD_KEY_H_
#define BAM_RECORD_KEY_H_

class WithBarcode {};
class NoBarcode {};

template<typename A, typename B>
using disable_if_same_or_derived =
typename std::enable_if<
    !std::is_base_of<A, typename
    std::remove_reference<B>::type
    >::value
>::type;

bool isRev(const seqan::BamAlignmentRecord &record)
{
    return (record.flag & 0x10) != 0;
}

template <typename THasBarcode = NoBarcode>
struct BamRecordKey
{
    BamRecordKey init(const seqan::BamAlignmentRecord &record) noexcept
    {
        pos = static_cast<uint64_t>(record.rID) << 32 |
            static_cast<uint64_t>((record.beginPos + static_cast<uint64_t>((isRev(record) == true ? length(record.seq) : 0)))) << 1 |
            static_cast<uint64_t>(isRev(record));
        return *this;
    }
    BamRecordKey init(const unsigned int rID, const unsigned int pos5end, const bool reverseStrand) noexcept
    {
        pos = static_cast<uint64_t>(rID) << 32 |
            static_cast<uint64_t>(pos5end) << 1 |
            static_cast<uint64_t>(reverseStrand);
        return *this;
    }
    BamRecordKey(const uint64_t pos) : pos(pos) {};

    template <typename TRecord, typename X =
        disable_if_same_or_derived<BamRecordKey, TRecord >>
    BamRecordKey(TRecord&& record)
    {
        init(std::forward<TRecord>(record));
    }

    __int32 get5EndPosition() const
    {
        return static_cast<__int32>(pos) >> 1;
    }

    __int32 getRID() const
    {
        return static_cast<__int32>(pos >> 32);
    }

    bool isReverseStrand() const
    {
        return (pos & 0x01) != 0;
    }
    bool friend operator<(const BamRecordKey<NoBarcode>& lhs, const BamRecordKey<NoBarcode>& rhs)
    {
        return lhs.pos < rhs.pos;
    }
    bool friend operator<=(const BamRecordKey<NoBarcode>& lhs, const BamRecordKey<NoBarcode>& rhs)
    {
        return lhs.pos <= rhs.pos;
    }
    bool friend lessEqualWithoutStrand(const BamRecordKey<THasBarcode>& rhs, const BamRecordKey<THasBarcode>& lhs)
    {
        return (rhs.pos & 0xFFFFFFFFFFFFFFFE) <= (lhs.pos & 0xFFFFFFFFFFFFFFFE);
    }
    bool friend operator==(const BamRecordKey<THasBarcode>& rhs, const BamRecordKey<THasBarcode>& lhs)
    {
        return rhs.pos == lhs.pos;
    }
    bool friend isEqualWithoutStrand(const BamRecordKey<THasBarcode>& rhs, const BamRecordKey<THasBarcode>& lhs)
    {
        return (rhs.pos & 0xFFFFFFFFFFFFFFFE) == (lhs.pos & 0xFFFFFFFFFFFFFFFE);
    }
private:
    uint64_t pos;
};

template <>
struct BamRecordKey<WithBarcode> : BamRecordKey<NoBarcode>
{
    template <typename TRecord>
    BamRecordKey(TRecord&& record) : BamRecordKey<NoBarcode>(record)
    {
        const std::string idString = toCString(record.qName);
        auto posStart = idString.find("TL:");
        if (posStart == std::string::npos)
        {
            barcode.clear();
            return;
        }
        posStart += 3;
        auto posEnd = idString.find(':', posStart);
        if (posEnd == std::string::npos)
            posEnd = idString.length();
        barcode = idString.substr(posStart, posEnd - posStart);
    }
    bool friend operator<(const BamRecordKey<WithBarcode>& lhs, const BamRecordKey<WithBarcode>& rhs)
    {
        if (operator<(static_cast<const BamRecordKey<NoBarcode>&>(lhs), static_cast<const BamRecordKey<NoBarcode>&>(rhs)))
            return true;
        if ((operator==(static_cast<const BamRecordKey<NoBarcode>&>(lhs), static_cast<const BamRecordKey<NoBarcode>&>(rhs))))
            return lhs.barcode < rhs.barcode;
        return false;
    }
private:
    std::string barcode;
};

// returns false if keys are from different chromosomes
template <typename THasBarcode>
bool calculate5EndDistance(const BamRecordKey<THasBarcode>& key1, const BamRecordKey<THasBarcode>& key2, int& distance)
{
    if (key1.getRID() != key2.getRID())
    {
        //distance = 0;
        return false;
    }
    distance = key2.get5EndPosition() - key1.get5EndPosition();
    return true;
}

#endif
