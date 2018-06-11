#include <deque>
#include <sstream>
#include <algorithm>

#include "dbg.h"


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::process_read(const string &dna_str, const uint32_t color_id, bool phase_first) {
    // the first KMER-sized substring will be $$$...$
    bitset<KMERBITS> akmer = 0, pkmer = 0;

    if (!phase_first && !second_phase_started) {
        do_sampling();
        second_phase_started = true;
    }

    // read the rest of the dna string...
    if (phase_first) {
        add_new_node(akmer, true, 255);
        for (uint64_t i = 0, sl = dna_str.length(); i < sl; ++i) {
            uint8_t sid = symbol_to_id(dna_str[i]);
            if (sid != 255) {
                akmer >>= LOGSIGMA;
                akmer |= shifted_sids[sid];
                add_new_node(akmer, false, sid);
            }
        }
        // the last edge leads to $
        dbg_kmers[akmer] |= 1;
    }
    else {
        // we will store the color for the outgoing edges of $$$...$
        colors[akmer];
        for (uint64_t i = 0, sl = dna_str.length(); i < sl; ++i) {
            uint8_t sid = symbol_to_id(dna_str[i]);
            if (sid != 255) {
                akmer >>= LOGSIGMA;
                akmer |= shifted_sids[sid];
                if (colors.find(pkmer) != colors.end()) {
                    colors[pkmer][sid].add(color_id);
                }
                pkmer = akmer;
            }
        }
        colors[pkmer][0].add(color_id);
    }
}


template<uint16_t KMERBITS>
inline void DeBrujinGraph<KMERBITS>::add_new_node(const bitset<KMERBITS> &akmer, bool new_node, uint8_t pc) {
    static bitset<KMERBITS> pkmer;
    if (!new_node) {
        dbg_kmers[pkmer] |= 1 << pc;
    }

    // Note that if there is no akmer in dbg_kmers, it will get allocated
    dbg_kmers[akmer];

    pkmer = akmer;
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::gen_succinct_dbg(const string& fname) {
    cerr << "Generating Succinct De Bruijn Graph..." << endl;
    cerr << "Sorting edge list..." << endl;

    // map<bitset<KMERBITS>, uint8_t, DeBrujinGraph<KMERBITS>::compare_lexicographically<KMERBITS>> dbg_kmers_sorted;

    typedef typename stxxl::VECTOR_GENERATOR<pair<bitset<KMERBITS>, uint8_t>>::result dbg_kmer_vector_type;
    dbg_kmer_vector_type dbg_kmers_sorted;
    dbg_kmers_sorted.resize(dbg_kmers.size());
    for (auto it = dbg_kmers.cbegin(), i = 0; it != dbg_kmers.end(); ++it, ++i) {
        dbg_kmers_sorted[i] = make_pair(it->first, it->second);
    }

    std::sort(dbg_kmers_sorted.begin(), dbg_kmers_sorted.end(),
            [this](const pair<bitset<KMERBITS>, uint8_t>& a, const pair<bitset<KMERBITS>, uint8_t>& b) -> bool {
                for (int i = kmer_bits - 1; i >= 0; --i) {
                    if (a.first[i] ^ b.first[i]) {
                        return b.first[i];
                    }
                }
                return false;
            });


    // for (auto it = dbg_kmers.cbegin(); it != dbg_kmers.end(); ++it) {
    //     dbg_kmers_sorted[it->first] = it->second;
    // }
    // dbg_kmers.clear();
    //
    cerr << "Saving edge list to file..." << endl;
    char buffer[LOGSIGMA * 8000] = { };
    size_t buffer_index = 0;
    ofstream f(fname, ios::out | ios::binary);
    auto save_data = [&f, &buffer, &buffer_index](size_t size){
        f.write(buffer, size);
        std::fill(buffer, buffer + sizeof(buffer), 0);
        buffer_index = 0;
    };
    for (size_t i = 0; i < dbg_kmers_sorted.size(); ++i) {
        for (uint8_t j = 0; j < SIGMA + 1; ++j) {
            if (((1 << j) & dbg_kmers_sorted[i].second) != 0) {
                uint8_t ibits = id_to_bits(j);
                for (uint8_t k = 0; k < LOGSIGMA; ++k, ++buffer_index) {
                    if ((1 << k) & ibits) {
                        buffer[buffer_index / 8] |= (1 << (buffer_index % 8));
                    }
                }
                if ((buffer_index / 8) >= sizeof(buffer)) {
                    save_data(sizeof(buffer));
                }
            }
        }
    }
    if (buffer_index > 0) {
        save_data(buffer_index);
    }

    auto gen_bin_list = [&dbg_kmers_sorted, &buffer, &buffer_index, &save_data](auto fn) {
        for (auto it = dbg_kmers_sorted.cbegin(); it != dbg_kmers_sorted.cend(); ++it) {
            // const auto &anode = nodes_ingoing[dbg_kmers.rank(it)];

            uint8_t ic = fn(it->first);
            if (ic > 0) {
                buffer_index += ic - 1;
                if (buffer_index >= sizeof(buffer)) {
                    size_t pbi = buffer_index;
                    save_data(sizeof(buffer));
                    buffer_index = pbi % sizeof(buffer);
                }
                buffer[buffer_index / 8] |= (1 << (buffer_index % 8));
            }
        }
        if (buffer_index > 0) {
            save_data(buffer_index);
        }
    };

    cerr << "Generating B_F list..." << endl;
    gen_bin_list([this](auto kmer) { return indegree(kmer); });

    cerr << endl << "Generating B_L list..." << endl;
    gen_bin_list([&dbg_kmers_sorted, this](auto kmer) { return outdegree(dbg_kmers[kmer]); });

    f.close();
}


template<typename KT, typename VT>
static size_t sparse_hash_map_difference(const sparse_hash_map<KT, VT> &left, const sparse_hash_map<KT, VT> &right) {
    size_t s = 0;
    for (auto it = left.cbegin(); it != left.cend(); ++it) {
        if (right.find(it->first) == right.end()) {
            ++s;
        }
    }

    return s;
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::do_sampling() {
    // maximum distance without
    uint32_t max_distance = max(sampling_max_distance, (uint32_t) log2(dbg_kmers.size()));
    cerr << "Starting sampling process with max distance: " << max_distance << "..." << std::endl;
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    for (auto it = dbg_kmers.cbegin(); it != dbg_kmers.cend(); ++it) {
        uint32_t outdeg = outdegree(dbg_kmers[it->first]);
        if (!visited[it->first] && (outdeg > 1)) { // || indegree(it->first) > 1   || it == dbg_kmers.cbegin()
            // store the colours on each branching node...
            colors[it->first];
            explicitly_stored_colors += outdeg;

            // if we do sampling (0 means that we only store colour information in branching nodes...)
            if (sampling_max_distance > 0) {
                for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                    if (((1 << i) & it->second) != 0) {
                        size_t path_length = 1;

                        bitset<KMERBITS> akmer = it->first;
                        akmer >>= LOGSIGMA;
                        akmer |= shifted_sids[i];
                        while (dbg_kmers.find(akmer) != dbg_kmers.end() && outdegree(dbg_kmers[akmer]) == 1) {
                            // && indegree(akmer) == 1
                            visited[akmer] = 1;
                            auto &anode = dbg_kmers[akmer];
                            for (uint8_t j = 0; j < SIGMA + 1; ++j) {
                                if (((1 << j) & anode) != 0) {
                                    akmer >>= LOGSIGMA;
                                    akmer |= shifted_sids[j];
                                    ++path_length;

                                    // store extra colour information (sampling)
                                    if (path_length % max_distance == 0) {
                                        colors[akmer];
                                        explicitly_stored_colors++;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        visited[it->first] = 1;
    }
}


template<uint16_t KMERBITS>
inline uint8_t DeBrujinGraph<KMERBITS>::outdegree(const bitset<SIGMA + 1> &ar) {
    uint8_t s = 0;
    for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        s += ar[i];
    }
    return s;
}


template<uint16_t KMERBITS>
inline uint8_t DeBrujinGraph<KMERBITS>::indegree(bitset<KMERBITS> pkmer) {
    static const bitset<KMERBITS> mask(string(kmer_bits, '1'));
    uint8_t s = 0;
    // uint64_t akmer = pkmer;
    uint8_t ac = bits_to_id((uint8_t) (pkmer >> (kmer_bits - LOGSIGMA)).to_ulong());
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


template<uint16_t KMERBITS>
inline bitset<MAXCOLORS> DeBrujinGraph<KMERBITS>::color_to_bitset(const Roaring &rc) {
    bitset<MAXCOLORS> acolor;
    for (uint32_t j = 0; j < C; ++j) {
        acolor[j] = rc.contains(j);
    }
    return acolor;
}


// // define static variables
// template<uint16_t KMERBITS>
// uint16_t DeBrujinGraph<KMERBITS>::kmer_bits = 0;


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


static inline char bits_to_char(const uint8_t s) {
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


static inline uint8_t bits_to_id(const uint8_t s) {
    switch (s) {
        case 0b000:
            return 0;
        case 0b001:
            return 1;
        case 0b011:
            return 2;
        case 0b101:
            return 3;
        case 0b111:
            return 4;
        default:
            return 255;
    }
}


static inline uint8_t id_to_bits(const uint8_t i) {
    switch (i) {
        case 0:
            return 0b000;
        case 1:
            return 0b001;
        case 2:
            return 0b011;
        case 3:
            return 0b101;
        case 4:
            return 0b111;
        default:
            return 255;
    }
}


static inline uint8_t symbol_to_id(const char c) {
    switch (c) {
        case '$':
            return 0;
        case 'A':
        case 'a':
            return 1;
        case 'C':
        case 'c':
            return 2;
        case 'G':
        case 'g':
            return 3;
        case 'T':
        case 't':
            return 4;
        default:
            return 255;
    }
}
