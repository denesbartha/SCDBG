#include <iostream>
#include <iomanip>
#include <string>
#include <sparsepp/spp.h>
#include <iostream>
#include <fstream>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/tuple/size.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>

#include "dbg.hpp"

using spp::sparse_hash_map;

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
    }
    else {
        std::cerr << "Unable to open file '" << fname << "'" << std::endl;
        exit(1);
    }
}


int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "usage: scdbg <kmer> <color> <input file> <output file>" << std::endl;
        std::cerr << "example usage: ./scdbg 32 100 input_file output_file" << std::endl;
        exit(1);
    }

    uint8_t kmer_size = (uint8_t) std::stoul(argv[1]);
    uint32_t colors = (uint32_t)std::stoi(argv[2]);

    // const uint8_t KMER = 3;
    // const uint8_t COLORS = 1;
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

    // const uint8_t KMER = 3;
    // const uint8_t COLORS = 1;
    // DeBrujinGraph<uint128_t, std::map<uint128_t, size_t>, std::map<uint128_t, size_t>::iterator> dbg(KMER);
    DeBrujinGraph<> dbg(3, 1);
    std::string dna_str1 = "TACGTCGACGACT";
    std::string dna_str2 = "TACGCGACT";
    // result:
    // TCCGTGGGACTAAA$C
    //  001111110111111    <-- B_F
    // 1110111100111111    <-- B_L
    dbg.process_read(dna_str1, 0);
    dbg.process_read(dna_str2, 1);
    dbg.process_read(dna_str1, 0, false);
    dbg.process_read(dna_str2, 1, false);
    dbg.do_stats();
    dbg.gen_succinct_dbg(argv[4]);


    // // DeBrujinGraph<uint128_t, std::set<uint128_t>, std::vector<uint128_t>::iterator> dbg(KMER);
    // // DeBrujinGraph<boost::multiprecision::uint128_t, std::map<boost::multiprecision::uint128_t, size_t>, std::map<boost::multiprecision::uint128_t, size_t>::iterator> dbg(kmer_size, colors);
    // // DeBrujinGraph<boost::multiprecision::uint256_t, std::map<boost::multiprecision::uint256_t, size_t>, std::map<boost::multiprecision::uint256_t, size_t>::iterator> dbg(kmer_size, colors);
    // DeBrujinGraph<uint64_t, std::map<uint64_t, size_t>, std::map<uint64_t, size_t>::iterator> dbg(kmer_size, colors);
    // for (int i = 0; i < 2; ++i) {
    //     std::cerr << "Start phase " << i << std::endl;
    //     std::ifstream myfile(argv[3]);
    //     std::string fname;
    //
    //     uint32_t color_id = 0;
    //     if (myfile.is_open()) {
    //         while (getline(myfile, fname)) {
    //             std::string dna_str = parse_inputfile(fname);
    //             // std::cout << dna_str.size() << std::endl;
    //             dbg.process_read(dna_str, color_id, i == 0);
    //             if (++color_id >= colors) {
    //                 break;
    //             }
    //         }
    //         myfile.close();
    //
    //         if (i != 0) {
    //             dbg.do_stats();
    //             // dbg.gen_succinct_dbg(argv[4]);
    //         }
    //     }
    //     else {
    //         std::cout << "Unable to open file '" << argv[1] << "'" << std::endl;
    //         exit(1);
    //     }
    //
    // }
    return 0;
}
 