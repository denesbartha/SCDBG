#include <iostream>
#include <iomanip>
#include <string>
#include <sparsepp/spp.h>
#include <iostream>
#include <fstream>

#include "dbg_wrapper.hpp"

using spp::sparse_hash_map;


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


    // DBGWrapper dbg(2, 1);
    // string dna_str = "TAATGCCATGGGATGTT";
    // dbg.process_read(dbg, dna_str, 0);
    // dbg.do_stats(dbg);


    // DBGWrapper dbg(3, 2);
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
    // dbg.gen_succinct_dbg(argv[5]);


    uint32_t colors = (uint32_t)stoi(argv[2]);
    uint32_t psmd = (uint32_t)stoi(argv[3]);
    DBGWrapper dbg(kmer_size, colors, psmd);
    for (int i = 0; i < 2; ++i) {
        cerr << "Start phase " << i << endl;
        ifstream myfile(argv[4]);
        string fname;

        uint32_t color_id = 0;
        if (myfile.is_open()) {
            cerr << "processing... " << endl;
            while (getline(myfile, fname)) {
                string dna_str = parse_inputfile(fname);
                // cout << dna_str.size() << endl;
                cerr << color_id << " ";
                dbg.process_read(dna_str, color_id, i == 0);
                if (++color_id >= colors) {
                    break;
                }
            }
            cerr << endl;
            myfile.close();

            if (i != 0) {
                // dbg.do_stats();
                dbg.gen_succinct_dbg(argv[5]);
            }
        }
        else {
            cout << "Unable to open file '" << argv[1] << "'" << endl;
            exit(1);
        }
    }
    return 0;
}
 