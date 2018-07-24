#ifndef SCDBG_UTILS_HPP
#define SCDBG_UTILS_HPP

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


inline void
buffer_to_file(std::ostream& f, char *data_buffer, size_t size, size_t& buffer_index, bool reset_index = true) {
    f.write(data_buffer, size);
    // reset the buffer
    std::fill(data_buffer, data_buffer + sizeof(data_buffer), 0);
    if (reset_index) {
        buffer_index = 0;
    }
}


inline size_t divide_and_to_upper(size_t a, size_t b) {
    return a % b == 0 ? a / b : (a / b) + 1;
}


inline void buffer_to_file(std::ostream& f, char *data_buffer, size_t length) {
    f.write(data_buffer, length);
}


inline void load_from_file(std::ifstream&f, char*& buffer, size_t length) {
    auto bytes = divide_and_to_upper(length, 8);
    buffer = new char[bytes];
    f.read(buffer, bytes);
}

inline void reset_buffer(char* buffer) {
    if (buffer != nullptr) {
        delete[] buffer;
        buffer = nullptr;
    }
}

#endif //SCDBG_UTILS_HPP
