#ifndef STATISTICS_DBG_H
#define STATISTICS_DBG_H

#include <set>
#include <array>
#include <bitset>
#include <sstream>
#include <algorithm>
#include <roaring.c>
#include <roaring.hh>
#include <sparsepp/spp.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

using spp::sparse_hash_map;

// SIGMA (for DNA, it is 4...)
#define SIGMA       4

//log_2(SIGMA + 1)
#define LOGSIGMA    3

// Node struct that contains nodes and the outgoing edges (this is necessary for the building process...)
template<uint8_t COLORBITS>
struct Node {
    std::bitset<SIGMA + 1> out_edges;
    // std::bitset<COLORBITS> colors;
};


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


template<uint8_t KMER, uint8_t COLORBITS, uint16_t KMERBITS = LOGSIGMA * KMER>
class DeBrujinGraph {
public:
    void process_read(const std::string &dna_str, const uint64_t color_id) {
        static uint8_t sid;
        std::cerr << "process " << color_id << "..." << std::endl;
        // read the first KMER-sized substring will be $$$..$
        std::bitset<KMERBITS> akmer = 0;

        // read the rest of the dna string...
        add_new_node(akmer, color_id, true, false, 255);
        for (uint64_t i = 0, sl = dna_str.length(); i < sl; ++i) {
            sid = symbol_to_bits(dna_str[i]);
            akmer >>= LOGSIGMA;
            akmer |= sid << (KMERBITS - LOGSIGMA);
            add_new_node(akmer, color_id, false, i == sl - 1, symbol_to_id(dna_str[i]));
        }
    }

    void do_stats() {
        std::cout << "Creating statistics..." << std::endl;

        uint64_t num_of_nodes = dbg.size();
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
        sparse_hash_map<std::bitset<COLORBITS>, std::bitset<1> > colors;

        for (const auto &item : dbg) {
            // uint64_t icnt = count_edges(item.second->in_edges);
            uint64_t icnt = count_incoming_edges(item.first);
            uint64_t ocnt = count_outgoing_edges(item.second.out_edges);
            in_degrees[icnt]++;
            out_degrees[ocnt]++;

            num_of_edges += icnt + ocnt;

            // colors[item.second.colors] = 1;

            print_node(item.first, item.second, icnt, ocnt);

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
        std::cout << "k-mer size:\t\t\t" << (uint64_t) KMER << std::endl;
        std::cout << "# of nodes:\t\t\t" << num_of_nodes << std::endl;
        std::cout << "# of edges:\t\t\t" << (num_of_edges / 2) << std::endl;
        std::cout << "# of colors:\t\t\t" << (uint64_t) COLORBITS << std::endl;
        std::cout << "# of color classes:\t\t" << colors.size() << std::endl;
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
        bool operator()(const std::bitset<N>& a, const std::bitset<N>& b) {
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
        std::cerr << "Sorting the edges..." << std::endl;

        std::set<std::bitset<KMERBITS + SIGMA + 1>, compare_lexicographically<KMERBITS + SIGMA + 1> > sdbg;
        std::bitset<KMERBITS + SIGMA + 1> node_edge;
        // for (auto it = dbg.cbegin(); it != dbg.cend(); ++it) {
        while (!dbg.empty()) {
            auto it = dbg.cbegin();
            if (KMERBITS <= 64) {
                node_edge = std::bitset<KMERBITS + SIGMA + 1>(it->first.to_ullong());
            }
            else {
                node_edge = std::bitset<KMERBITS + SIGMA + 1>(it->first.to_string());
            }
            node_edge <<= SIGMA + 1;
            node_edge |= it->second.out_edges.to_ulong();
            sdbg.insert(node_edge);
            dbg.erase(it);
        }


        std::cerr << "Generating output..." << std::endl;
        std::bitset<LOGSIGMA> edge_label;
        std::bitset<KMERBITS> node;
        // for (auto it = sdbg.cbegin(); it != sdbg.cend(); ++it) {
        while (!sdbg.empty()) {
            auto it = sdbg.cbegin();
            std::bitset<SIGMA + 1> edges = (uint8_t)it->to_ullong();
            if (KMERBITS <= 64) {
                node = std::bitset<KMERBITS>((*it >> (SIGMA + 1)).to_ullong());
            }
            else {
                node = std::bitset<KMERBITS>((*it >> (SIGMA + 1)).to_string());
            }

            // std::cout << it->to_string() << " " << kmer_to_str(node) << " " << edges << std::endl;

            uint8_t edges_bits = (uint8_t)it->to_ulong();
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (((1 << i) & edges_bits) != 0) {
                    std::cout << base[i];
                }
            }
            // std::cout << std::endl;

            sdbg.erase(it);
        }
        std::cout << std::endl;
    }


private:
    inline void add_new_node(const std::bitset<KMERBITS> &akmer, const uint64_t color_id, bool new_node, bool end_node, uint8_t pc) {
        static std::bitset<KMERBITS> pkmer;
        if (!new_node) {
            dbg[pkmer].out_edges.set(pc);
        }

        // Note that if there is no akmer in dbg, it will be allocated
        Node<COLORBITS> *anode = &(dbg[akmer]);
        // add_color(anode->colors, color_id);
        if (end_node) {
            anode->out_edges[0] = 1;
        }
        pkmer = akmer;
    }

    inline uint8_t count_outgoing_edges(const std::bitset<SIGMA + 1> &ar) {
        uint8_t s = 0;
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            s += ar[i];
        }
        return s;
    }

    // TODO: this implementation works only for SIGMA=4...
    inline uint8_t count_incoming_edges(const std::bitset<KMERBITS> &pkmer) {
        uint8_t s = 0;
        std::bitset<KMERBITS> akmer = pkmer;
        uint8_t ac = symbol_to_id(bits_to_char((uint8_t) (akmer >> (KMERBITS - LOGSIGMA)).to_ulong()));
        akmer <<= LOGSIGMA;
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            akmer ^= symbol_to_bits(base[i]);
            if (dbg.contains(akmer) && dbg[akmer].out_edges[ac]) {
                ++s;
            }
            akmer ^= symbol_to_bits(base[i]);
        }
        return s;
    }


    // from a given bitstring generates the appropriate kmer string
    std::string kmer_to_str(const std::bitset<KMERBITS> &kmer_str) {
        std::stringstream ss;
        // ss << str.to_string() << std::endl;
        for (int i = 0; i < KMER; ++i) {
            uint8_t ac = 0;
            for (int j = LOGSIGMA - 1; j >= 0; --j) {
                ac <<= 1;
                ac |= kmer_str[i * LOGSIGMA + j];
            }
            ss << bits_to_char(ac);
        }
        return ss.str();
    }


    void print_node(const std::bitset<KMERBITS> &str, const Node<COLORBITS> &node, uint64_t icnt, uint64_t ocnt) {
        std::cout << kmer_to_str(str) << std::endl;

        // print the in edges
        std::cout << std::endl << "in edges:  ";
        std::bitset<KMERBITS> akmer = str;
        uint8_t ac = symbol_to_id(bits_to_char((uint8_t) (akmer >> (KMERBITS - LOGSIGMA)).to_ulong()));
        akmer <<= LOGSIGMA;
        for (int i = 0; i < SIGMA + 1; ++i) {
            akmer ^= symbol_to_bits(base[i]);
            if (dbg.contains(akmer) && dbg[akmer].out_edges[ac]) {
                std::cout << base[i] << " ";
            }
            // std::cout << base[i] << ": " << ((dbg.contains(akmer)) ? dbg[akmer].out_edges[ac] : 0) << std::endl;
            akmer ^= symbol_to_bits(base[i]);
        }


        // print the out edges
        std::cout << std::endl << "out edges: ";
        for (int i = 0; i < SIGMA + 1; ++i) {
            if (node.out_edges[i]) {
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

    static inline void add_color(std::bitset<COLORBITS> &colors, uint64_t acolor) {
        colors.set(acolor);
    }

    sparse_hash_map<std::bitset<KMERBITS>, Node<COLORBITS>> dbg;
    // number of edges
    uint64_t M = 0;
    const char base[5] = {'$', 'A', 'C', 'G', 'T'};

    // typedef int vt[10];
    // constexpr const vt val() { return {0, 1, 2, 3}; }
    // const int vv[] = val();

    // Roaring64Map r2;
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
