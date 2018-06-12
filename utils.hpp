#ifndef SCDBG_UTILS_H
#define SCDBG_UTILS_H

#include <iostream>
#include <map>
#include <sparsepp/spp.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

template<typename A, typename B>
std::pair<B, A> flip_pair(const std::pair<A, B> &p) {
    return std::pair<B, A>(p.second, p.first);
}


template<typename A, typename B>
std::multimap<B, A> flip_map(const spp::sparse_hash_map<A, B> &src) {
    std::multimap<B, A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), flip_pair<A, B>);
    return dst;
}


inline uint8_t log2_int(size_t x) {
    uint8_t i = 0;
    do {
        ++i;
    } while ((x >>= 1) != 0);
    return i;
}

#endif //SCDBG_UTILS_H
