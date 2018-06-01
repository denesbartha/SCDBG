#ifndef SCDBG_BITSET_WRAPPER_H
#define SCDBG_BITSET_WRAPPER_H

#include <set>
#include <climits>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

// typedef __int128 int128_t;
// typedef unsigned __int128 uint128_t;

template<typename T, typename BitSetClass, typename Iterator_Type>
class bitset_wrapper {
public:
    inline size_t size() const;

    inline void add(const T);

    inline size_t rank(const T);

    inline bool contains(const T) const;

    Iterator_Type begin() { return dbg_kmers.begin(); }

    Iterator_Type end() { return dbg_kmers.end(); }

private:
    BitSetClass dbg_kmers;

    bool first_rank = true;
};

template<>
uint64_t bitset_wrapper<uint64_t, Roaring64Map, Roaring64MapSetBitForwardIterator>::size() const {
    return dbg_kmers.cardinality();
}

template<>
size_t bitset_wrapper<uint64_t, std::map<uint64_t, size_t>, std::map<uint64_t, size_t>::iterator>::size() const {
    return dbg_kmers.size();
}

template<>
size_t bitset_wrapper<boost::multiprecision::uint128_t, std::map<boost::multiprecision::uint128_t, size_t>,
        std::map<boost::multiprecision::uint128_t, size_t>::iterator>::size() const {
    return dbg_kmers.size();
}

template<>
size_t
bitset_wrapper<boost::multiprecision::uint256_t, std::map<boost::multiprecision::uint256_t, size_t>,
        std::map<boost::multiprecision::uint256_t, size_t>::iterator>::size() const {
    return dbg_kmers.size();
}


template<>
void bitset_wrapper<uint64_t, Roaring64Map, Roaring64MapSetBitForwardIterator>::add(const uint64_t x) {
    dbg_kmers.add(x);
}

template<>
void bitset_wrapper<uint64_t, std::map<uint64_t, size_t>, std::map<uint64_t, size_t>::iterator>::add(
        const uint64_t x) {
    dbg_kmers[x] = 0;
}

template<>
void
bitset_wrapper<boost::multiprecision::uint128_t, std::map<boost::multiprecision::uint128_t, size_t>,
        std::map<boost::multiprecision::uint128_t, size_t>::iterator>::add(
        const boost::multiprecision::uint128_t x) {
    dbg_kmers[x] = 0;
}

template<>
void bitset_wrapper<boost::multiprecision::uint256_t, std::map<boost::multiprecision::uint256_t, size_t>,
        std::map<boost::multiprecision::uint256_t, size_t>::iterator>::add(const boost::multiprecision::uint256_t x) {
    dbg_kmers[x] = 0;
}


template<>
uint64_t bitset_wrapper<uint64_t, Roaring64Map, Roaring64MapSetBitForwardIterator>::rank(const uint64_t x) {
    return dbg_kmers.rank(x) - 1;
}

template<>
size_t
bitset_wrapper<uint64_t, std::map<uint64_t, size_t>, std::map<uint64_t, size_t>::iterator>::rank(const uint64_t x) {
    // for the first time calling the rank operation - calculate the indexes of the entries of dbg_kmers
    if (first_rank) {
        size_t i = 0;
        for (auto it = dbg_kmers.begin(); it != dbg_kmers.end(); ++it) {
            it->second = i++;
        }
        first_rank = false;
    }
    return dbg_kmers[x];
}

template<>
size_t bitset_wrapper<boost::multiprecision::uint128_t, std::map<boost::multiprecision::uint128_t, size_t>,
        std::map<boost::multiprecision::uint128_t, size_t>::iterator>::rank(const boost::multiprecision::uint128_t x) {
    // for the first time calling the rank operation - calculate the indexes of the entries of dbg_kmers
    if (first_rank) {
        size_t i = 0;
        for (auto it = dbg_kmers.begin(); it != dbg_kmers.end(); ++it) {
            it->second = i++;
        }
        first_rank = false;
    }
    return dbg_kmers[x];
}

template<>
size_t bitset_wrapper<boost::multiprecision::uint256_t, std::map<boost::multiprecision::uint256_t, size_t>,
        std::map<boost::multiprecision::uint256_t, size_t>::iterator>::rank(const boost::multiprecision::uint256_t x) {
    // for the first time calling the rank operation - calculate the indexes of the entries of dbg_kmers
    if (first_rank) {
        size_t i = 0;
        for (auto it = dbg_kmers.begin(); it != dbg_kmers.end(); ++it) {
            it->second = i++;
        }
        first_rank = false;
    }
    return dbg_kmers[x];
}


template<>
bool bitset_wrapper<uint64_t, Roaring64Map, Roaring64MapSetBitForwardIterator>::contains(const uint64_t x) const {
    return dbg_kmers.contains(x);
}

template<>
bool bitset_wrapper<uint64_t, std::map<uint64_t, size_t>, std::map<uint64_t, size_t>::iterator>::contains(
        const uint64_t x) const {
    return dbg_kmers.find(x) != dbg_kmers.end();
}

template<>
bool
bitset_wrapper<boost::multiprecision::uint128_t, std::map<boost::multiprecision::uint128_t, size_t>,
        std::map<boost::multiprecision::uint128_t, size_t>::iterator>::contains(
        const boost::multiprecision::uint128_t x) const {
    return dbg_kmers.find(x) != dbg_kmers.end();
}

template<>
bool bitset_wrapper<boost::multiprecision::uint256_t, std::map<boost::multiprecision::uint256_t, size_t>,
        std::map<boost::multiprecision::uint256_t, size_t>::iterator>::contains(
        const boost::multiprecision::uint256_t x) const {
    return dbg_kmers.find(x) != dbg_kmers.end();
}

#endif //SCDBG_BITSET_WRAPPER_H
