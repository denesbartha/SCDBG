#ifndef STATISTICS_DBG_H
#define STATISTICS_DBG_H

#include <set>
#include <array>
#include <bitset>
#include <climits>
#include <sstream>
#include <algorithm>
#include <roaring.c>
#include <roaring.hh>
#include <stxxl/vector>
#include <sparsepp/spp.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

using namespace std;
using spp::sparse_hash_map;

// SIGMA (for DNA, it is 4...)
#define SIGMA       4

//log_2(SIGMA + 1)
#define LOGSIGMA    3

#define MAXDISTANCE 100

static inline uint8_t symbol_to_bits(const char c);
static inline char bits_to_char(uint8_t s);
static inline uint8_t symbol_to_id(const char c);
static const char base[5] = {'$', 'A', 'C', 'G', 'T'};


template<uint8_t KMERBITS>
class DeBrujinGraph {
public:

    DeBrujinGraph() : DeBrujinGraph(32, 10, 100) { }

    DeBrujinGraph(const uint8_t pkm, const uint32_t pc, const uint32_t psmd) : km(pkm),
                                                                               C(pc), sampling_max_distance(psmd) {
        kmer_bits = LOGSIGMA * pkm;

        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            bitset<KMERBITS> sid = symbol_to_bits(base[i]);
            sid <<= kmer_bits - LOGSIGMA;
            shifted_sids[i] = sid;
        }
    }

    void process_read(const string &dna_str, const uint32_t color_id, bool phase_first = true);

    void do_stats();


    void gen_succinct_dbg(string fname);


private:
    inline void add_new_node(const bitset<KMERBITS> &akmer, const uint32_t color_id, bool new_node, uint8_t pc);

    inline uint8_t outdegree(const bitset<SIGMA + 1> &ar);


    inline uint8_t indegree(bitset<KMERBITS> pkmer);

    string kmer_to_str(bitset<KMERBITS> kmer_str);

    void print_node(const bitset<KMERBITS> &str, uint64_t icnt, uint64_t ocnt);

    void do_sampling();

    array<bitset<KMERBITS>, SIGMA + 1> shifted_sids;

    uint8_t km;
    static uint16_t kmer_bits;
    uint32_t sampling_max_distance;
    size_t explicitly_stored_colors = 0;

    // functor for comparing two given bitsets lexicographically

    template<size_t N>
    struct compare_lexicographically {
        bool operator()(const bitset<N> &x, const bitset<N> &y) {
            for (int i = kmer_bits - 1; i >= 0; --i) {
                if (x[i] ^ y[i]) {
                    return y[i];
                }
            }
            return false;
        }
    };

    // sparse_hash_map<bitset<KMERBITS>, uint8_t> dbg_kmers;
    sparse_hash_map<bitset<KMERBITS>, array<Roaring, SIGMA + 1>> colors;
    map<bitset<KMERBITS>, uint8_t, compare_lexicographically<KMERBITS>> dbg_kmers;

    // number of edges
    size_t M = 0;
    // number of colors
    uint32_t C;

    bool second_phase_started = false;

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
