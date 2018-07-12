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

    DeBrujinGraph() : DeBrujinGraph(32, 10, -1) {}

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

    void sort_dbg();

    void save_edge_list(ofstream& f);

    void save_table_F(ofstream& f);

    void save_colors(const string& fname);

    void save_color_classes(const string& fname, const multimap<size_t, bitset<MAXCOLORS>>& ordered_cm,
                            sparse_hash_map<bitset<MAXCOLORS>, size_t>& color_class_order, size_t set_bits);

    void save_store_vector(ostream& f);

    void save_color_bit_vector(ostream& f, sparse_hash_map<bitset<MAXCOLORS>, size_t>& color_class_order,
                               bool boundary);

protected:

    inline uint8_t outdegree(const bitset<SIGMA + 1>& ar);

    inline uint8_t indegree(bitset<KMERBITS> pkmer);

    inline bitset<MAXCOLORS> color_to_bitset(const Roaring& rc);

    array<bitset<KMERBITS>, SIGMA + 1> shifted_sids;

    uint8_t km;
    uint32_t kmer_bits;
    uint32_t sampling_max_distance;
    size_t explicitly_stored_colors = 0;

    sparse_hash_map<bitset<KMERBITS>, uint8_t> dbg_kmers;
    sparse_hash_map<bitset<KMERBITS>, array<Roaring, SIGMA + 1>> colors;
    // map<bitset<KMERBITS>, uint8_t, compare_lexicographically<KMERBITS>> dbg_kmers;

    // number of edges
    size_t num_of_edges = 0;
    // number of colors
    uint32_t num_of_colors;

private:
    bool second_phase_started = false;

    typedef typename stxxl::VECTOR_GENERATOR<pair<bitset<KMERBITS>, uint8_t>>::result dbg_kmer_vector_type;
    dbg_kmer_vector_type dbg_kmers_sorted;

    // Roaring64Map dbg_kmers;
    // bitset_wrapper<T, BitSetClass, Iterator_Type> dbg_kmers;

    // typedef typename stxxl::VECTOR_GENERATOR<bitset<SIGMA + 1>>::result vector;
    // typedef typename vector<bitset<SIGMA + 1>> vector;
    // vector nodes_outgoing;

    // typedef typename stxxl::VECTOR_GENERATOR<array<Roaring, SIGMA + 1>>::result color_vector_type;
    // color_vector_type colors;

    // vector<array<Roaring, SIGMA + 1>> colors;
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
