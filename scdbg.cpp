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


string parse_inputfile(const string &fname) {
    ifstream myfile(fname);
    string line;
    if (myfile.is_open()) {
        // the first line is just a header...
        getline(myfile, line);

        // the second line contains the actual DNA data (flatten...)
        getline(myfile, line);

        myfile.close();
        return line;
    }
    else {
        cerr << "Unable to open file '" << fname << "'" << endl;
        exit(1);
    }
}


class DBGWrapper {
public:
    DBGWrapper(const uint8_t pkmer_size, const uint32_t pcolors, const uint32_t psmd) : kmer_size(pkmer_size) {
        if (kmer_size <= KMER8BYTES) {
            dbg8 = new DeBrujinGraph<64>(pkmer_size, pcolors, psmd);
        }
        else if (kmer_size <= KMER16BYTES) {
            dbg16 = new DeBrujinGraph<128>(pkmer_size, pcolors, psmd);
        }
        else if (kmer_size <= KMER24BYTES) {
            dbg24 = new DeBrujinGraph<192>(pkmer_size, pcolors, psmd);
        }
        else if (kmer_size <= KMER32BYTES) {
            dbg32 = new DeBrujinGraph<256>(pkmer_size, pcolors, psmd);
        }
        else if (kmer_size <= KMER40BYTES) {
            dbg40 = new DeBrujinGraph<320>(pkmer_size, pcolors, psmd);
        }
        else {
            cerr << "Maximal k-mer size is " << KMER40BYTES << "..." << endl;
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

    void process_read(const string &dna_str, const uint32_t color_id, bool phase_first) {
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

    void gen_succinct_dbg(string fname) {
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
};


int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "usage: scdbg <kmer> <color> <sampling distance> <input file> <output file>" << endl;
        cerr << "example usage: ./scdbg 32 100 100 input_file output_file" << endl;
        exit(1);
    }

    uint8_t kmer_size = (uint8_t) stoul(argv[1]);
    if (kmer_size > KMER40BYTES) {
        cerr << "Maximal k-mer size is " << KMER40BYTES << "..." << endl;
        exit(1);
    }

    uint32_t colors = (uint32_t)stoi(argv[2]);
    uint32_t psmd = (uint32_t)stoi(argv[3]);


    // DBGWrapper dbg(kmer_size, colors);
    // string dna_str = "TAATGCCATGGGATGTT";
    // dbg.process_read(dbg, dna_str, 0);
    // dbg.do_stats(dbg);


    // DBGWrapper dbg(3, 1, 100);
    // string dna_str1 = "TACGTCGACGACT";
    // string dna_str2 = "TACGCGACT";
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


    DBGWrapper dbg(kmer_size, colors, psmd);
    for (int i = 0; i < 2; ++i) {
        cerr << "Start phase " << i << endl;
        ifstream myfile(argv[4]);
        string fname;

        uint32_t color_id = 0;
        if (myfile.is_open()) {
            while (getline(myfile, fname)) {
                string dna_str = parse_inputfile(fname);
                // cout << dna_str.size() << endl;
                dbg.process_read(dna_str, color_id, i == 0);
                if (++color_id >= colors) {
                    break;
                }
            }
            myfile.close();

            if (i != 0) {
                dbg.do_stats();
                // dbg.gen_succinct_dbg(argv[5]);
            }
        }
        else {
            cout << "Unable to open file '" << argv[1] << "'" << endl;
            exit(1);
        }

    }
    return 0;
}
 