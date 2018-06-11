#ifndef SCDBG_DBG_WRAPPER_H
#define SCDBG_DBG_WRAPPER_H

#include "dbg.h"
#include "dbg.cpp"
#include "dbg_stats.h"
#include "dbg_stats.cpp"

#define KMER8BYTES 64 / LOGSIGMA
#define KMER16BYTES 128 / LOGSIGMA
#define KMER24BYTES 192 / LOGSIGMA
#define KMER32BYTES 256 / LOGSIGMA
#define KMER40BYTES 320 / LOGSIGMA


class DBGWrapper {
public:
    DBGWrapper(const uint8_t pkmer_size, const uint32_t pcolors, const uint32_t psmd = 0) : kmer_size(pkmer_size) {
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
            ((DeBrujinGraphStats<64>*)dbg8)->do_stats();
        }
        else if (kmer_size <= KMER16BYTES) {
            ((DeBrujinGraphStats<128>*)dbg16)->do_stats();
        }
        else if (kmer_size <= KMER24BYTES) {
            ((DeBrujinGraphStats<192>*)dbg24)->do_stats();
        }
        else if (kmer_size <= KMER32BYTES) {
            ((DeBrujinGraphStats<256>*)dbg32)->do_stats();
        }
        else if (kmer_size <= KMER40BYTES) {
            ((DeBrujinGraphStats<320>*)dbg40)->do_stats();
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

#endif //SCDBG_DBG_WRAPPER_H
