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

#include "bitset_wrapper.h"


using spp::sparse_hash_map;

// SIGMA (for DNA, it is 4...)
#define SIGMA       4

//log_2(SIGMA + 1)
#define LOGSIGMA    3


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


static inline char bits_to_char(uint8_t s) {
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


static inline uint8_t symbol_to_id(const char c) {
    switch (c) {
        case '$':
            return 0;
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 3;
        case 'T':
            return 4;
        default:
            return 255;
    }
}


template<typename T = uint64_t, typename BitSetClass = Roaring64Map,
        typename Iterator_Type = Roaring64MapSetBitForwardIterator>
class DeBrujinGraph {
public:

    DeBrujinGraph() : DeBrujinGraph(32) {}

    DeBrujinGraph(const uint8_t pkm) : km(pkm), kmer_bits(LOGSIGMA * pkm) {}

    void process_read(const std::string &dna_str, const uint64_t color_id, bool phase_first = true) {
        static T sid;
        std::cerr << "process " << color_id << "..." << std::endl;
        // read the first KMER-sized substring will be $$$..$
        // std::bitset<KMERBITS> akmer = 0;
        // boost::multiprecision::uint256_t akmer = 0;
        T akmer = 0;

        if (!phase_first && !second_phase_started) {
            nodes_outgoing.resize(dbg_kmers.size());
            for (uint64_t i = 0; i < nodes_outgoing.size(); ++i) {
                nodes_outgoing[i] = 0;
            }

            second_phase_started = true;
        }

        // read the rest of the dna string...
        if (phase_first) {
            dbg_kmers.add(akmer);
        }
        else {
            add_new_node(akmer, color_id, true, false, 255);
        }
        for (uint64_t i = 0, sl = dna_str.length(); i < sl; ++i) {
            sid = symbol_to_bits(dna_str[i]);
            akmer >>= LOGSIGMA;
            akmer |= sid << (kmer_bits - LOGSIGMA);

            if (sid != 255) {
                if (phase_first) {
                    dbg_kmers.add(akmer);
                }
                else {
                    add_new_node(akmer, color_id, false, i == sl - 1, symbol_to_id(dna_str[i]));
                }
            }
        }
    }

    void do_stats() {
        std::cout << "Creating statistics..." << std::endl;

        uint64_t num_of_nodes = dbg_kmers.size();
        uint64_t num_of_edges = 0;
        std::vector<uint64_t> in_degrees(SIGMA + 2, 0);
        std::vector<uint64_t> out_degrees(SIGMA + 2, 0);
        uint64_t in_out_one = 0;
        uint64_t branching_in = 0;
        uint64_t branching_out = 0;
        uint64_t branching_in_and_out = 0;
        uint64_t branching_in_or_out = 0;
        uint64_t source = 0;
        uint64_t sink = 0;
        // uint64_t
        // uint64_t colored_branches = 0;
        // sparse_hash_map<std::bitset<COLORBITS>, std::bitset<1> > colors;

        std::cerr << dbg_kmers.size() << std::endl;

        for (const auto &item : dbg_kmers) {
            // std::bitset<kmer_bits> bb = item;
            // std::cout << bb.to_string() << " " << nodes_outgoing[dbg_kmers.rank(item) - 1] << std::endl;

            // uint64_t icnt = count_edges(item.second->in_edges);
            uint64_t icnt = count_incoming_edges(item);
            uint64_t ocnt = count_outgoing_edges(nodes_outgoing[dbg_kmers.rank(item)]);
            in_degrees[icnt]++;
            out_degrees[ocnt]++;

            num_of_edges += icnt + ocnt;

            // colors[item.second.colors] = 1;

            // print_node(item, icnt, ocnt);

            if (icnt == 1 && ocnt == 1) {
                in_out_one++;
            }
            else if (icnt > 1 && ocnt <= 1) {
                branching_in++;
                branching_in_or_out++;
            }
            else if (ocnt > 1 && icnt <= 1) {
                branching_out++;
                branching_in_or_out++;
            }
            else if (icnt > 1 && ocnt > 1) {
                branching_in_and_out++;
                branching_in_or_out++;
            }

            // it can be source/sink and branching node at the same time...
            if (icnt == 0) {
                source++;
            }
            if (ocnt == 0) {
                sink++;
            }

            // delete item.second;
        }

        std::cout << "in degrees:" << std::endl;
        for (int i = 0; i < SIGMA + 1; ++i) {
            std::cout << i << ": " << in_degrees[i] << std::endl;
        }

        std::cout << "out degrees:" << std::endl;
        for (int i = 0; i < SIGMA + 1; ++i) {
            std::cout << i << ": " << out_degrees[i] << std::endl;
        }

        std::cout << std::fixed;
        std::cout << std::setprecision(5);
        std::cout << "k-mer size:\t\t\t" << (uint64_t) km << std::endl;
        std::cout << "# of nodes:\t\t\t" << num_of_nodes << std::endl;
        std::cout << "# of edges:\t\t\t" << (num_of_edges / 2) << std::endl;
        // std::cout << "# of colors:\t\t\t" << (uint64_t) COLORBITS << std::endl;
        // std::cout << "# of color classes:\t\t" << colors.size() << std::endl;
        std::cout << "#source nodes:\t\t\t" << source << "\t\t\t" << (float) source / num_of_nodes << "%" << std::endl;
        std::cout << "#sink nodes:\t\t\t" << sink << "\t\t\t" << (float) sink / num_of_nodes << "%" << std::endl;
        std::cout << "#in, out = 1:\t\t\t" << in_out_one << "\t\t\t" << (float) in_out_one / num_of_nodes << "%"
                  << std::endl;
        std::cout << "#brancing in > 1:\t\t" << branching_in << "\t\t\t" << (float) branching_in / num_of_nodes << "%"
                  << std::endl;
        std::cout << "#brancing out > 1:\t\t" << branching_out << "\t\t\t" << (float) branching_out / num_of_nodes
                  << "%"
                  << std::endl;
        std::cout << "#brancing in and out > 1:\t" << branching_in_and_out << "\t\t\t"
                  << (float) branching_in_and_out / num_of_nodes << "%" << std::endl;
        std::cout << "#brancing in or out > 1:\t" << branching_in_or_out << "\t\t\t"
                  << (float) branching_in_or_out / num_of_nodes << "%" << std::endl;
    }


    // functor for comparing two given bitsets lexicographically
    template<std::size_t N>
    struct compare_lexicographically {
        bool operator()(const std::bitset<N> &a, const std::bitset<N> &b) {
            for (int i = N - 1; i >= 0; i--) {
                if (a[i] ^ b[i]) {
                    return b[i];
                }
            }
            return false;
        }
    };


    void gen_succinct_dbg() {
        std::cerr << "Generating Succinct De Bruijn Graph..." << std::endl;
        std::cout << "Generating Edge list..." << std::endl;
        for (auto it : dbg_kmers) {
            const auto &anode = nodes_outgoing[dbg_kmers.rank(it)];
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (anode[i]) {
                    std::cout << base[i];
                }
            }
        }
        std::cout << std::endl;

        std::cout << "Generating B_F list..." << std::endl;
        for (auto it : dbg_kmers) {
            // const auto &anode = nodes_ingoing[dbg_kmers.rank(it)];
            int ic = count_incoming_edges(it);
            if (ic > 0) {
                std::cout << std::string(ic - 1, '0') << "1";
            }
            else {
                std::cout << " ";
            }
        }

        std::cout << std::endl << "Generating B_L list..." << std::endl;
        for (auto it : dbg_kmers) {
            const auto &anode = nodes_outgoing[dbg_kmers.rank(it)];
            std::cout << std::string(count_outgoing_edges(anode) - 1, '0') << "1";
        }
    }


private:
    inline void add_new_node(const T akmer, const uint64_t color_id, bool new_node, bool end_node, uint8_t pc) {
        // static std::bitset<KMERBITS> pkmer;
        static uint64_t pkmer_rank;
        if (!new_node) {
            // dbg[pkmer].out_edges.set(pc);
            nodes_outgoing[pkmer_rank].set(pc);
        }

        // Note that if there is no akmer in dbg, it will be allocated
        // add_color(anode->colors, color_id);

        pkmer_rank = dbg_kmers.rank(akmer);
        auto &anode = nodes_outgoing[pkmer_rank];
        // if (!new_node) {
        //     nodes_incoming[pkmer_rank].set();
        // }
        if (end_node) {
            // anode->out_edges[0] = 1;
            anode.set(0);
        }
    }

    inline uint8_t count_outgoing_edges(const std::bitset<SIGMA + 1> &ar) {
        uint8_t s = 0;
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            s += ar[i];
        }
        return s;
    }


    inline uint8_t count_incoming_edges(T pkmer) {
        static T mask = (((T)1 << kmer_bits) - 1);
        uint8_t s = 0;
        // uint64_t akmer = pkmer;
        uint8_t ac = symbol_to_id(bits_to_char((uint8_t) (pkmer >> (kmer_bits - LOGSIGMA))));
        pkmer = (pkmer << LOGSIGMA) & mask;
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            pkmer ^= symbol_to_bits(base[i]);
            if (dbg_kmers.contains(pkmer) && nodes_outgoing[dbg_kmers.rank(pkmer)][ac]) {
                ++s;
            }
            pkmer ^= symbol_to_bits(base[i]);
        }
        return s;
    }


    // from a given bitstring generates the appropriate kmer string
    std::string kmer_to_str(const T kmer_str) {
        std::stringstream ss;
        // ss << str.to_string() << std::endl;
        for (int i = 0; i < km; ++i) {
            uint8_t ac = 0;
            for (int j = LOGSIGMA - 1; j >= 0; --j) {
                ac <<= 1;
                ac |= ((kmer_str & ((T)1 << (i * LOGSIGMA + j))) != 0) ? 1 : 0;
            }
            ss << bits_to_char(ac);
        }
        return ss.str();
    }


    void print_node(const T str, uint64_t icnt, uint64_t ocnt) {
        static T mask = (((T)1 << kmer_bits) - 1);
        std::cout << std::bitset<31 * LOGSIGMA>(str).to_string() << std::endl << kmer_to_str(str) << std::endl;

        // print the in edges
        std::cout << std::endl << "in edges:  ";
        T akmer = str;
        // the last added character
        uint8_t ac = symbol_to_id(bits_to_char((uint8_t) (akmer >> (kmer_bits - LOGSIGMA))));
        // generate all the possible (SIGMA + 1) k-mers that could be the source
        akmer = (akmer << LOGSIGMA) & mask;
        for (int i = 0; i < SIGMA + 1; ++i) {
            akmer ^= symbol_to_bits(base[i]);
            if (dbg_kmers.contains(akmer) && nodes_outgoing[dbg_kmers.rank(akmer)][ac]) {
                std::cout << base[i] << " ";
            }
            akmer ^= symbol_to_bits(base[i]);
        }


        // print the out edges
        std::cout << std::endl << "out edges: ";
        const auto &anode = nodes_outgoing[dbg_kmers.rank(str)];
        for (int i = 0; i < SIGMA + 1; ++i) {
            if (anode[i]) {
                std::cout << base[i] << " ";
            }
        }
        std::cout << std::endl;

        if (icnt == 1 && ocnt == 1) {
            std::cout << "in out, 1" << std::endl;
        }
        else if (icnt > 1 && ocnt <= 1) {
            std::cout << "branching in > 1" << std::endl;
        }
        else if (ocnt > 1 && icnt <= 1) {
            std::cout << "branching out > 1" << std::endl;
        }
        else if (icnt > 1 && ocnt > 1) {
            std::cout << "branching in and out > 1" << std::endl;
        }

        if (icnt == 0) {
            std::cout << "source" << std::endl;
        }
        else if (ocnt == 0) {
            std::cout << "sink" << std::endl;
        }

        std::cout << std::endl;
    }

    // static inline void add_color(std::bitset<COLORBITS> &colors, uint64_t acolor) {
    //     colors.set(acolor);
    // }

    uint8_t km;
    uint16_t kmer_bits;

    // sparse_hash_map<std::bitset<KMERBITS>, Node<COLORBITS>> dbg;


    // Roaring64Map dbg_kmers;
    bitset_wrapper<T, BitSetClass, Iterator_Type> dbg_kmers;

    // typedef typename stxxl::VECTOR_GENERATOR<std::bitset<SIGMA + 1>>::result vector;
    typedef typename std::vector<std::bitset<SIGMA + 1>> vector;
    vector nodes_outgoing;
    // vector nodes_ingoing;

    // number of edges
    size_t M = 0;
    const char base[5] = {'$', 'A', 'C', 'G', 'T'};

    // typedef int vt[10];
    // constexpr const vt val() { return {0, 1, 2, 3}; }
    // const int vv[] = val();

    bool second_phase_started = false;

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
