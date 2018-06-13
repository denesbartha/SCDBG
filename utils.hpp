#ifndef SCDBG_UTILS_H
#define SCDBG_UTILS_H

#include <ostream>
#include <map>
#include <sparsepp/spp.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

static inline uint8_t symbol_to_bits(const char c) {
    switch (c) {
        case '$':
            return 0b000;
        case 'A':
            return 0b001;
        case 'C':
            return 0b011;
        case 'G':
            return 0b101;
        case 'T':
            return 0b111;
        default:
            return 255;
    }
}


static inline char bits_to_char(const uint8_t s) {
    switch (s) {
        case 0b000:
            return '$';
        case 0b001:
            return 'A';
        case 0b011:
            return 'C';
        case 0b101:
            return 'G';
        case 0b111:
            return 'T';
        default:
            return 0;
    }
}


static inline uint8_t bits_to_id(const uint8_t s) {
    switch (s) {
        case 0b000:
            return 0;
        case 0b001:
            return 1;
        case 0b011:
            return 2;
        case 0b101:
            return 3;
        case 0b111:
            return 4;
        default:
            return 255;
    }
}


static inline uint8_t id_to_bits(const uint8_t i) {
    switch (i) {
        case 0:
            return 0b000;
        case 1:
            return 0b001;
        case 2:
            return 0b011;
        case 3:
            return 0b101;
        case 4:
            return 0b111;
        default:
            return 255;
    }
}


static inline uint8_t symbol_to_id(const char c) {
    switch (c) {
        case '$':
            return 0;
        case 'A':
        case 'a':
            return 1;
        case 'C':
        case 'c':
            return 2;
        case 'G':
        case 'g':
            return 3;
        case 'T':
        case 't':
            return 4;
        default:
            return 255;
    }
}


template<typename A, typename B>
std::pair<B, A> flip_pair(const std::pair<A, B>& p) {
    return std::pair<B, A>(p.second, p.first);
}


template<typename A, typename B>
std::multimap<B, A> flip_map(const spp::sparse_hash_map<A, B>& src) {
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


inline void
buffer_to_file(std::ostream& f, char *data_buffer, size_t size, size_t& buffer_index, bool reset_index = true) {
    f.write(data_buffer, size);
    // reset the buffer
    std::fill(data_buffer, data_buffer + sizeof(data_buffer), 0);
    if (reset_index) {
        buffer_index = 0;
    }
}


inline void buffer_to_file(std::ostream& f, char *data_buffer, size_t size) {
    f.write(data_buffer, size);
}

#endif //SCDBG_UTILS_H
