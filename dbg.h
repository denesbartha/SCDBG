#ifndef STATISTICS_DBG_H
#define STATISTICS_DBG_H

#include <set>
#include <array>
#include <bitset>
#include <roaring.c>
#include <roaring.hh>
#include <stxxl/vector>
#include <sparsepp/spp.h>

#include "utils.hpp"
#include "config.h"

using namespace std;
using spp::sparse_hash_map;


template<uint16_t KMERBITS>
class DeBrujinGraph {
public:

    DeBrujinGraph() : DeBrujinGraph(32, 10, 0) {}

    DeBrujinGraph(const uint8_t pkm, const uint32_t pc, const uint32_t psmd = 0) : km(pkm), kmer_bits(LOGSIGMA * pkm),
                                                                                   sampling_max_distance(psmd),
                                                                                   num_of_colors(pc) {
        // kmer_bits = (uint16_t)LOGSIGMA * pkm;

        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            bitset<KMERBITS> sid = symbol_to_bits(base[i]);
            sid <<= kmer_bits - LOGSIGMA;
            shifted_sids[i] = sid;
        }
    }

    void process_read(const string& dna_str, const uint32_t color_id, bool phase_first);

    void gen_succinct_dbg(const string& fname);

private:
    inline void add_new_node(const bitset<KMERBITS>& akmer, bool new_node, uint8_t pc);

    void do_sampling();

    void save_edge_list(ofstream& f);

    void save_table_F(ofstream& f);

    void save_colors(const string& fname);

    void save_color_classes(const string& fname, const multimap<size_t, bitset<MAXCOLORS>>& ordered_cm,
                            sparse_hash_map<bitset<MAXCOLORS>, pair<size_t, uint8_t>>& color_class_order,
                            size_t set_bits);

    void save_store_vector(ostream& f);

    void save_color_bit_vector(const string& fname, sparse_hash_map<bitset<MAXCOLORS>, pair<size_t, uint8_t>>& color_class_order,
                               bool boundary);

protected:

    inline uint8_t outdegree(const bitset<SIGMA + 1>& ar);

    inline uint8_t indegree(bitset<KMERBITS> pkmer);

    // inline bitset<MAXCOLORS> color_to_bitset(const Roaring& rc);

    inline void add_color(size_t& kmer_color_hash, const uint32_t color_id);

    inline size_t add_color_class(const bitset<MAXCOLORS>& bitvector);

    void sort_dbg();

    struct compare_bit_vector : public std::less<pair<bitset<KMERBITS>, uint8_t>> {
        const pair<bitset<KMERBITS>, uint8_t> bs;

        const pair<bitset<KMERBITS>, uint8_t> mask = pair<bitset<KMERBITS>, uint8_t>(bitset<KMERBITS>(string(LOGSIGMA, '1')), 255);

        int kmer_bits;

        compare_bit_vector(int pkmer_bits) : kmer_bits(pkmer_bits) { }

        bool operator()(const pair<bitset<KMERBITS>, uint8_t>& a, const pair<bitset<KMERBITS>, uint8_t>& b) const {
            for (int i = kmer_bits - 1; i >= 0; --i) {
                if (a.first[i] ^ b.first[i]) {
                    return b.first[i];
                }
            }
            return false;
        }

        pair<bitset<KMERBITS>, uint8_t> min_value() const { return bs; }

        pair<bitset<KMERBITS>, uint8_t> max_value() const { return mask; }
    };

    array<bitset<KMERBITS>, SIGMA + 1> shifted_sids;

    uint8_t km;
    uint32_t kmer_bits;
    uint32_t sampling_max_distance;
    size_t explicitly_stored_colors = 0;

    sparse_hash_map<bitset<KMERBITS>, uint8_t> dbg_kmers;

    sparse_hash_map<bitset<KMERBITS>, array<size_t, SIGMA + 1>> colors;

    struct color_class_t {
        color_class_t() {}

        color_class_t(const bitset<MAXCOLORS>& pbitvector) : bitvector(pbitvector) {}

        bitset<MAXCOLORS> bitvector;
        size_t cnt = 0;
    };

    sparse_hash_map<uint64_t, color_class_t> color_classes;
    hash<bitset<MAXCOLORS>> hash_color_class;
    hash<size_t> hash_int;
    size_t set_bits = 0;

    size_t num_of_edges = 0;
    uint32_t num_of_colors;

protected:
    bool second_phase_started = false;

    typedef typename stxxl::VECTOR_GENERATOR<pair<bitset<KMERBITS>, uint8_t>>::result dbg_kmer_vector_type;
    dbg_kmer_vector_type dbg_kmers_sorted;
    // vector<pair<bitset<KMERBITS>, uint8_t>> dbg_kmers_sorted;

};


constexpr uint8_t calc_edge_cnt(const uint8_t edge) {
    uint8_t s = 0;
    for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        if (((1 << i) & edge) != 0) {
            ++s;
        }
    }
    return s;
}


#endif //STATISTICS_DBG_H
