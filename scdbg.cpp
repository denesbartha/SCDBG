#include <iostream>
#include <iomanip>
#include <string>
#include <sparsepp/spp.h>
#include <iostream>
#include <fstream>

#define MAXCOLOR        32
#define ITERATIONSTART  65
#define ITERATIONSEND   66  // 32 * 32

#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/tuple/size.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include "dbg.hpp"

using spp::sparse_hash_map;



#define PRED(r, state) BOOST_PP_NOT_EQUAL( \
    BOOST_PP_TUPLE_ELEM(2, 0, state) \
  , BOOST_PP_TUPLE_ELEM(2, 1, state) \
)

#define OP(r, state) ( \
    BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)) \
  , BOOST_PP_TUPLE_ELEM(2, 1, state) \
)

// calculates the upper_bound(log_2(c))
constexpr uint8_t calc_colorbits(uint64_t c) {
    uint8_t s = 0;
    if (c == 0) {
        return 0;
    }

    for (uint8_t i = 1; i <= 64; ++i) {
        if ((c & 1) == 1) {
            s = i;
        }
        c >>= 1;
    }
    return s;
}

// kmer length
// #define KMER    10

// the number of colors
// #define COLORS  95146

// log(2 * SIGMA)
#define LOG_2SIGMA  3

// #define KMERBITS    (KMER * 2)

// SIGMA #numbers
// typedef uint64_t sarray[SIGMA];

std::string parse_inputfile(const std::string &fname) {
    std::ifstream myfile(fname);
    std::string line;
    if (myfile.is_open()) {
        // the first line is just a header...
        getline(myfile, line);

        // the second line contains the actual DNA data (flatten...)
        getline(myfile, line);

        myfile.close();
        return line;
    } else {
        std::cerr << "Unable to open file '" << fname << "'" << std::endl;
        exit(1);
    }
}


#define MACRO(r, state) \
    case BOOST_PP_TUPLE_ELEM(2, 0, state): { \
        DeBrujinGraph<BOOST_PP_TUPLE_ELEM(2, 0, state) / MAXCOLOR, \
                      (BOOST_PP_TUPLE_ELEM(2, 0, state) % MAXCOLOR)> dbg; \
        std::string dna_str = "TAATGCCATGGGATGTT"; \
        dbg.process_read(dna_str, 0); \
        dbg.do_stats(); \
        dbg.gen_succinct_dbg(); \
        break; }


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "usage: stats <input file>" << std::endl;
        exit(1);
    }

    const uint8_t KMER = 3;
    const uint8_t COLORS = 1;

    // switch (KMER * MAXCOLOR + COLORS) {
    //     BOOST_PP_FOR((ITERATIONSTART, ITERATIONSEND), PRED, OP, MACRO)
    //
    //     default:
    //         std::cerr << "Invalid k-mer or colour..." << std::endl;
    // }
    //
    // DeBrujinGraph<KMER, COLORS> dbg;
    // std::string dna_str = "TAATGCCATGGGATGTT";
    // process_read(dbg, dna_str, 0);
    // do_stats(dbg);

    DeBrujinGraph<KMER, COLORS> dbg;
    std::string dna_str1 = "TACGTCGACGACT";
    std::string dna_str2 = "TACGCGACT";
    // result:
    // TCCGTGGGACTAAA$C
    //  001111110111111    <-- B_F
    // 1110111100111111    <-- B_L
    dbg.process_read(dna_str1, 0);
    dbg.process_read(dna_str2, 1);
    dbg.do_stats();
    dbg.gen_succinct_dbg();



    // DeBrujinGraph<KMER, COLORS> dbg;
    // std::ifstream myfile(argv[1]);
    // std::string fname;
    // uint64_t color_id = 0;
    // if (myfile.is_open()) {
    //     while (getline(myfile, fname)) {
    //         std::string dna_str = parse_inputfile(fname);
    //         dbg.process_read(dna_str, color_id);
    //         if (++color_id >= COLORS) {
    //             break;
    //         }
    //     }
    //
    //     myfile.close();
    //
    //     dbg.do_stats();
    // } else {
    //     std::cout << "Unable to open file '" << argv[1] << "'" << std::endl;
    //     exit(1);
    // }

    return 0;
}
 