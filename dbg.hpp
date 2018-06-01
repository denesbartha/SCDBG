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


template <uint8_t KMERBITS>
class DeBrujinGraph {
public:

    DeBrujinGraph() : DeBrujinGraph(32, 10) { }

    std::array<std::bitset<KMERBITS>, SIGMA + 1> shifted_sids;
    DeBrujinGraph(const uint8_t pkm, const uint32_t pc) : km(pkm), kmer_bits(LOGSIGMA * pkm), C(pc) {
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            std::bitset<KMERBITS> sid = symbol_to_bits(base[i]);
            sid <<= kmer_bits - LOGSIGMA;
            shifted_sids[i] = sid;
        }
    }

    void process_read(const std::string &dna_str, const uint32_t color_id, bool phase_first = true) {
        // static T sid;
        std::cerr << "process " << color_id << "..." << std::endl;
        // read the first KMER-sized substring will be $$$..$
        // boost::multiprecision::uint256_t akmer = 0;
        std::bitset<KMERBITS> akmer = 0;

        if (!phase_first && !second_phase_started) {
            // nodes_outgoing.resize(dbg_kmers.size(), 0);
            // colors.resize(dbg_kmers.size());
            // for (uint64_t i = 0; i < dbg_kmers.size(); ++i) {
            //     nodes_outgoing[i] = 0;
            //
            // }

            second_phase_started = true;
        }

        // read the rest of the dna string...
        if (phase_first) {
            add_new_node(akmer, color_id, true, 255);
            for (uint64_t i = 0, sl = dna_str.length(); i < sl; ++i) {
                uint8_t sid = symbol_to_id(dna_str[i]);
                if (sid != 255) {
                    akmer >>= LOGSIGMA;
                    akmer |= shifted_sids[sid];
                    add_new_node(akmer, color_id, false, symbol_to_id(dna_str[i]));
                }
            }
            // the last edge leads to $...
            dbg_kmers[akmer] |= 1;
        }
    }

    void do_stats() {
        std::cerr << "Creating statistics..." << std::endl;

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

        for (const auto &item : dbg_kmers) {
            // std::bitset<kmer_bits> bb = item;
            // std::cout << bb.to_string() << " " << nodes_outgoing[dbg_kmers.rank(item) - 1] << std::endl;

            // uint64_t icnt = count_edges(item.second->in_edges);
            uint64_t icnt = count_incoming_edges(item.first);
            uint64_t ocnt = count_outgoing_edges(dbg_kmers[item.first]);
            in_degrees[icnt]++;
            out_degrees[ocnt]++;

            num_of_edges += icnt + ocnt;

            // colors[item.second.colors] = 1;

            // print_node(item.first, icnt, ocnt);

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
                print_node(item.first, icnt, ocnt);
            }
            if (ocnt == 0) {
                sink++;
            }

            // delete item.second;
        }

        std::cerr << "Generate color stats..." << std::endl;
        std::map<std::bitset<KMERBITS>, size_t> cm;

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
        std::cout << "# of colors:\t\t\t" << (uint64_t) C << std::endl;
        std::cout << "# of color classes:\t\t" << cm.size() << std::endl;
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


    void gen_succinct_dbg(std::string fname) {
        std::cerr << "Generating Succinct De Bruijn Graph..." << std::endl;
        std::cerr << "Generating Edge list..." << std::endl;
        for (auto it : dbg_kmers) {
            // std::cout << it.first.to_string() << std::endl;
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (((1 << i) & it.second) != 0) {
                    std::cout << base[i];
                }
            }
        }

        // std::cerr << "Generating B_F list..." << std::endl;
        // for (auto it : dbg_kmers) {
        //     // const auto &anode = nodes_ingoing[dbg_kmers.rank(it)];
        //     int ic = count_incoming_edges(it);
        //     if (ic > 0) {
        //         std::cout << std::string(ic - 1, '0') << "1";
        //     }
        //     else {
        //         std::cout << " ";
        //     }
        // }
        //
        // std::cerr << std::endl << "Generating B_L list..." << std::endl;
        // for (auto it : dbg_kmers) {
        //     const auto &anode = nodes_outgoing[dbg_kmers.rank(it)];
        //     std::cout << std::string(count_outgoing_edges(anode) - 1, '0') << "1";
        // }
    }


private:
    inline void add_new_node(const std::bitset<KMERBITS>& akmer, const uint32_t color_id, bool new_node, uint8_t pc) {
        static std::bitset<KMERBITS> pkmer;
        if (!new_node) {
            // colors[pkmer_rank][pc].add(color_id);
            dbg_kmers[pkmer] |= 1 << pc;
        }

        // Note that if there is no akmer in dbg_kmers, it will get allocated
        auto &anode = dbg_kmers[akmer];
        // if (end_node) {
        //     // anode->out_edges[0] = 1;
        //     // anode.set(0);
        //
        //     // colors[pkmer_rank][pc].add(color_id);
        //     anode = 1;
        // }
        pkmer = akmer;
    }

    inline uint8_t count_outgoing_edges(const std::bitset<SIGMA + 1> &ar) {
        uint8_t s = 0;
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            s += ar[i];
        }
        return s;
    }


    inline uint8_t count_incoming_edges(std::bitset<KMERBITS> pkmer) {
        static const std::bitset<KMERBITS> mask(std::string(kmer_bits, '1'));
        uint8_t s = 0;
        // uint64_t akmer = pkmer;
        uint8_t ac = symbol_to_id(bits_to_char((uint8_t)(pkmer >> (kmer_bits - LOGSIGMA)).to_ulong()));
        pkmer = (pkmer << LOGSIGMA) & mask;
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            pkmer ^= symbol_to_bits(base[i]);
            // if (dbg_kmers.contains(pkmer) && nodes_outgoing[dbg_kmers.rank(pkmer)][ac]) {
            if (dbg_kmers.find(pkmer) != dbg_kmers.end() && (dbg_kmers[pkmer] & (1 << ac)) != 0) {
                ++s;
            }
            pkmer ^= symbol_to_bits(base[i]);
        }
        return s;
    }


    // from a given bitstring generates the appropriate kmer string
    std::string kmer_to_str(std::bitset<KMERBITS> kmer_str) {
        static const std::bitset<KMERBITS> mask(std::string(LOGSIGMA, '1'));
        std::stringstream ss;
        // ss << str.to_string() << std::endl;
        for (int i = 0; i < km; ++i) {
            uint8_t ac = (uint8_t)(kmer_str & mask).to_ulong();
            ss << bits_to_char(ac);
            kmer_str >>= LOGSIGMA;
        }
        return ss.str();
    }


    void print_node(const std::bitset<KMERBITS>& str, uint64_t icnt, uint64_t ocnt) {
        static const std::bitset<KMERBITS> mask(std::string(kmer_bits, '1'));
        std::cout << std::bitset<KMERBITS>(str).to_string() << std::endl << kmer_to_str(str) << std::endl;

        // print the in edges
        std::cout << std::endl << "in edges:  ";
        std::bitset<KMERBITS> akmer = str;
        // the last added character
        uint8_t ac = symbol_to_id(bits_to_char((uint8_t) (akmer >> (kmer_bits - LOGSIGMA)).to_ulong()));
        // generate all the possible (SIGMA + 1) k-mers that could be the source
        akmer = (akmer << LOGSIGMA) & mask;
        for (int i = 0; i < SIGMA + 1; ++i) {
            akmer ^= symbol_to_bits(base[i]);
            // if (dbg_kmers.contains(akmer) && nodes_outgoing[dbg_kmers.rank(akmer)][ac]) {
            if (dbg_kmers.find(akmer) != dbg_kmers.end() && (dbg_kmers[akmer] & (1 << ac)) != 0) {
                std::cout << base[i] << " ";
            }
            akmer ^= symbol_to_bits(base[i]);
        }

        // print the out edges
        std::cout << std::endl << "out edges: ";
        const auto &anode = dbg_kmers[str];
        for (int i = 0; i < SIGMA + 1; ++i) {
            if (((1 << i) & anode) != 0) {
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

    // functor for comparing two given bitsets lexicographically

    template<std::size_t N>
    struct compare_lexicographically {
        bool operator()(const std::bitset<N> &a, const std::bitset<N> &b) {
            for (int i = N - 1; i >= 0; --i) {
                if (a[i] ^ b[i]) {
                    return b[i];
                }
            }
            return false;
        }
    };
    sparse_hash_map<std::bitset<KMERBITS>, uint8_t> dbg_kmers;
    // std::map<std::bitset<KMERBITS>, uint8_t, compare_lexicographically<KMERBITS>> dbg_kmers;


    // Roaring64Map dbg_kmers;
    // bitset_wrapper<T, BitSetClass, Iterator_Type> dbg_kmers;

    // typedef typename stxxl::VECTOR_GENERATOR<std::bitset<SIGMA + 1>>::result vector;
    // typedef typename std::vector<std::bitset<SIGMA + 1>> vector;
    // vector nodes_outgoing;

    typedef typename stxxl::VECTOR_GENERATOR<std::array<Roaring, SIGMA + 1>>::result color_vector_type;
    color_vector_type colors;
    // std::vector<std::array<Roaring, SIGMA + 1>> colors;

    // number of edges
    size_t M = 0;
    uint32_t C;
    const char base[5] = {'$', 'A', 'C', 'G', 'T'};

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
