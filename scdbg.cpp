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

#define KMER8BYTES 64 / LOGSIGMA
#define KMER16BYTES 128 / LOGSIGMA
#define KMER24BYTES 192 / LOGSIGMA
#define KMER32BYTES 256 / LOGSIGMA
#define KMER40BYTES 320 / LOGSIGMA


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


class DBGWrapper {
public:
    DBGWrapper(uint8_t pkmer_size, uint32_t pcolors) : kmer_size(pkmer_size), colors(pcolors) {
        if (kmer_size <= KMER8BYTES) {
            dbg8 = new DeBrujinGraph<64>(kmer_size, colors);
        }
        else if (kmer_size <= KMER16BYTES) {
            dbg16 = new DeBrujinGraph<128>(kmer_size, colors);
        }
        else if (kmer_size <= KMER24BYTES) {
            dbg24 = new DeBrujinGraph<192>(kmer_size, colors);
        }
        else if (kmer_size <= KMER32BYTES) {
            dbg32 = new DeBrujinGraph<256>(kmer_size, colors);
        }
        else if (kmer_size <= KMER40BYTES) {
            dbg40 = new DeBrujinGraph<320>(kmer_size, colors);
        }
        else {
            std::cerr << "Maximal k-mer size is " << KMER40BYTES << "..." << std::endl;
        }
    }

    ~DBGWrapper() {
        if (dbg8 != nullptr) {
            delete dbg8;
        }
        if (dbg16 != nullptr) {
            delete dbg16;
        }
        if (dbg24 != nullptr) {
            delete dbg24;
        }
        if (dbg32 != nullptr) {
            delete dbg32;
        }
        if (dbg40 != nullptr) {
            delete dbg40;
        }
    }

    void process_read(const std::string &dna_str, const uint32_t color_id, bool phase_first) {
        if (kmer_size <= KMER8BYTES) {
            dbg8->process_read(dna_str, color_id, phase_first);
        }
        else if (kmer_size <= KMER16BYTES) {
            dbg16->process_read(dna_str, color_id, phase_first);
        }
        else if (kmer_size <= KMER24BYTES) {
            dbg24->process_read(dna_str, color_id, phase_first);
        }
        else if (kmer_size <= KMER32BYTES) {
            dbg32->process_read(dna_str, color_id, phase_first);
        }
        else if (kmer_size <= KMER40BYTES) {
            dbg40->process_read(dna_str, color_id, phase_first);
        }
    }

    void do_stats() {
        if (kmer_size <= KMER8BYTES) {
            dbg8->do_stats();
        }
        else if (kmer_size <= KMER16BYTES) {
            dbg16->do_stats();
        }
        else if (kmer_size <= KMER24BYTES) {
            dbg24->do_stats();
        }
        else if (kmer_size <= KMER32BYTES) {
            dbg32->do_stats();
        }
        else if (kmer_size <= KMER40BYTES) {
            dbg40->do_stats();
        }
    }

    void gen_succinct_dbg(std::string fname) {
        if (kmer_size <= KMER8BYTES) {
            dbg8->gen_succinct_dbg(fname);
        }
        else if (kmer_size <= KMER16BYTES) {
            dbg16->gen_succinct_dbg(fname);
        }
        else if (kmer_size <= KMER24BYTES) {
            dbg24->gen_succinct_dbg(fname);
        }
        else if (kmer_size <= KMER32BYTES) {
            dbg32->gen_succinct_dbg(fname);
        }
        else if (kmer_size <= KMER40BYTES) {
            dbg40->gen_succinct_dbg(fname);
        }
    }


private:
    DeBrujinGraph<64>* dbg8 = nullptr;
    DeBrujinGraph<128>* dbg16 = nullptr;
    DeBrujinGraph<192>* dbg24 = nullptr;
    DeBrujinGraph<256>* dbg32 = nullptr;
    DeBrujinGraph<320>* dbg40 = nullptr;

    uint8_t kmer_size;
    uint32_t colors;
};


int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "usage: scdbg <kmer> <color> <input file> <output file>" << std::endl;
        std::cerr << "example usage: ./scdbg 32 100 input_file output_file" << std::endl;
        exit(1);
    }

    uint8_t kmer_size = (uint8_t) std::stoul(argv[1]);
    if (kmer_size > KMER40BYTES) {
        std::cerr << "Maximal k-mer size is " << KMER40BYTES << "..." << std::endl;
        exit(1);
    }

    uint32_t colors = (uint32_t)std::stoi(argv[2]);


    // DBGWrapper dbg(kmer_size, colors);
    // std::string dna_str = "TAATGCCATGGGATGTT";
    // dbg.process_read(dbg, dna_str, 0);
    // dbg.do_stats(dbg);


    // DBGWrapper dbg(3, 1);
    // std::string dna_str1 = "TACGTCGACGACT";
    // std::string dna_str2 = "TACGCGACT";
    // // result:
    // // TCCGTGGGACTAAA$C
    // //  001111110111111    <-- B_F
    // // 1110111100111111    <-- B_L
    // dbg.process_read(dna_str1, 0, true);
    // dbg.process_read(dna_str2, 1, true);
    // dbg.process_read(dna_str1, 0, false);
    // dbg.process_read(dna_str2, 1, false);
    // dbg.do_stats();
    // dbg.gen_succinct_dbg(argv[4]);


    DBGWrapper dbg(kmer_size, colors);
    for (int i = 0; i < 2; ++i) {
        std::cerr << "Start phase " << i << std::endl;
        std::ifstream myfile(argv[3]);
        std::string fname;

        uint32_t color_id = 0;
        if (myfile.is_open()) {
            while (getline(myfile, fname)) {
                std::string dna_str = parse_inputfile(fname);
                // std::cout << dna_str.size() << std::endl;
                dbg.process_read(dna_str, color_id, i == 0);
                if (++color_id >= colors) {
                    break;
                }
            }
            myfile.close();

            if (i != 0) {
                dbg.do_stats();
                // dbg.gen_succinct_dbg(argv[4]);
            }
        }
        else {
            std::cout << "Unable to open file '" << argv[1] << "'" << std::endl;
            exit(1);
        }

    }
    return 0;
}
 