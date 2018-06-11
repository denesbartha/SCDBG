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
void DeBrujinGraph<KMERBITS>::do_stats() {
    cerr << "Creating statistics..." << endl;

    uint64_t num_of_nodes = dbg_kmers.size();
    uint64_t num_of_edges = 0;
    vector<uint64_t> in_degrees(SIGMA + 2, 0);
    vector<uint64_t> out_degrees(SIGMA + 2, 0);
    uint64_t in_out_one = 0;
    uint64_t branching_in = 0;
    uint64_t branching_out = 0;
    uint64_t branching_in_and_out = 0;
    uint64_t branching_in_or_out = 0;
    uint64_t source = 0;
    uint64_t sink = 0;
    uint64_t maxpath_length = 0;
    uint64_t branching_color_is_the_same[SIGMA + 2] = {};
    uint64_t branching_outdegrees[SIGMA + 2] = {};
    // uint64_t color_is_the_same_as_incoming = 0;
    // uint64_t branching_outdegrees_sum = 0;

    sparse_hash_map<bitset<MAXCOLORS>, size_t> cm;
    sparse_hash_map<bitset<MAXCOLORS>, size_t> cm_rest;
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    sparse_hash_map<size_t, size_t> paths;
    for (auto it = dbg_kmers.cbegin(); it != dbg_kmers.cend(); ++it) {
        uint64_t icnt = indegree(it->first);
        uint64_t ocnt = outdegree(dbg_kmers[it->first]);
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
            // print_node(it->first, icnt, ocnt);
        }
        if (ocnt == 0) {
            sink++;
        }


        // traverse the paths
        if ((visited.find(it->first) == visited.end() && ocnt > 1) || it == dbg_kmers.cbegin()) {
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
                                break;
                            }
                        }
                    }

                    if (paths.find(path_length) != paths.end()) {
                        paths[path_length]++;
                    }
                    else {
                        paths[path_length] = 1;
                    }
                    if (path_length > maxpath_length) {
                        maxpath_length = path_length;
                    }
                }
            }
        }
        visited[it->first] = 1;

        // color stats...
        // cerr << "Generate color stats..." << endl;
        if (colors.find(it->first) != colors.end()) {
            if (ocnt > 1) {
                // auto current_node_color = get_color(it->first);
                Roaring prev_color;
                bool first_node = true;
                bool color_is_the_same = true;
                for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                    if (((1 << i) & it->second) != 0) {
                        auto ac = colors[it->first][i];
                        cm[color_to_bitset(ac)]++;
                        if (!first_node) {
                            color_is_the_same &= prev_color == ac;
                        }
                        prev_color = ac;
                        first_node = false;
                        // if (ac == current_node_color) {
                        //     color_is_the_same_as_incoming++;
                        // }
                        // branching_outdegrees_sum++;
                    }
                }
                if (color_is_the_same) {
                    branching_color_is_the_same[ocnt]++;
                }
                branching_outdegrees[ocnt]++;
            }
            else if (icnt > 1) {
                // cm_rest[color_to_bitset(get_color(it->first))]++;
                for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                    if (((1 << i) & it->second) != 0) {
                        auto ac = colors[it->first][i];
                        cm_rest[color_to_bitset(ac)]++;
                    }
                }
            }
        }
        // else if (icnt > 1 || ocnt > 1) {
        //     cm_rest[color_to_bitset(get_color(it->first))]++;
        //
        //     if (ocnt > 1) {
        //
        //     }
        // }
    }


    cout << "in degrees:" << endl;
    for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        cout << i << ": " << in_degrees[i] << endl;
    }

    cout << "out degrees:" << endl;
    for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        cout << i << ": " << out_degrees[i] << endl;
    }

    cout << fixed;
    cout << setprecision(5);
    cout << "k-mer size:\t\t\t" << (uint64_t) km << endl;
    cout << "# of nodes:\t\t\t" << num_of_nodes << endl;
    cout << "# of edges:\t\t\t" << (num_of_edges / 2) << endl;
    cout << "# of colors:\t\t\t" << (uint64_t) C << endl;
    cout << "# of color classes:\t\t" << cm.size() << endl;
    cout << "# of color classes - not stored:\t" << sparse_hash_map_difference<bitset<MAXCOLORS>, size_t>(cm_rest, cm)
         << endl;
    for (uint8_t i = 2; i <= SIGMA + 1; ++i) {
        cout << "# of B_{*," << i << "} nodes, where the color is the same:\t" << branching_color_is_the_same[i] << "/"
             << branching_outdegrees[i] << endl;
    }
    // cout << "# of B_{*,+} nodes, where the outgoing color is the same as the node's color: "
    //      << color_is_the_same_as_incoming << "/" << branching_outdegrees_sum << endl;
    cout << "# of explicitly stored colors:\t" << explicitly_stored_colors << endl;
    cout << "#source nodes:\t\t\t" << source << "\t\t\t" << (float) source / num_of_nodes << "%" << endl;
    cout << "#sink nodes:\t\t\t" << sink << "\t\t\t" << (float) sink / num_of_nodes << "%" << endl;
    cout << "#in, out = 1:\t\t\t" << in_out_one << "\t\t\t" << (float) in_out_one / num_of_nodes << "%"
         << endl;
    cout << "#brancing in > 1:\t\t" << branching_in << "\t\t\t" << (float) branching_in / num_of_nodes << "%"
         << endl;
    cout << "#brancing out > 1:\t\t" << branching_out << "\t\t\t" << (float) branching_out / num_of_nodes
         << "%"
         << endl;
    cout << "#brancing in and out > 1:\t" << branching_in_and_out << "\t\t\t"
         << (float) branching_in_and_out / num_of_nodes << "%" << endl;
    cout << "#brancing in or out > 1:\t" << branching_in_or_out << "\t\t\t"
         << (float) branching_in_or_out / num_of_nodes << "%" << endl;

    cout << "maximal path legth:\t\t" << maxpath_length << " " << paths[maxpath_length] << endl;

    size_t ps = 0;
    size_t os = 0;
    for (auto it = paths.cbegin(); it != paths.cend(); ++it) {
        ps += it->first * it->second;
        os += it->second;
    }
    cout << "average path legth:\t\t" << (double) ps / os << endl;
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


// Returns the colour of a given node
template<uint16_t KMERBITS>
inline Roaring DeBrujinGraph<KMERBITS>::get_color(const bitset<KMERBITS> &pkmer) {
    static const bitset<KMERBITS> mask(string(kmer_bits, '1'));
    Roaring rcolor;
    deque<bitset<KMERBITS>> kmer_queue(1, pkmer);
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    while (!kmer_queue.empty()) {
        auto akmer = kmer_queue.front();
        kmer_queue.pop_front();
        std::vector<bitset<KMERBITS>> incoming_nodes;
        // step backwards until we find the color information or hit a branching node
        uint8_t ac = 255;
        while (colors.find(akmer) == colors.end() && indegree(akmer) <= 1) {
            ac = bits_to_id((uint8_t) (akmer >> (kmer_bits - LOGSIGMA)).to_ulong());
            // cerr << kmer_to_str(akmer) << " " << (akmer >> (kmer_bits - LOGSIGMA)).to_ulong() << endl;
            akmer = (akmer << LOGSIGMA) & mask;
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                akmer ^= symbol_to_bits(base[i]);
                if (dbg_kmers.find(akmer) != dbg_kmers.end() && (dbg_kmers[akmer] & (1 << ac)) != 0) {
                    // there is at most only one incoming node => break if we have found that
                    break;
                }
                akmer ^= symbol_to_bits(base[i]);
            }
        }
        // if we have found the color => merge it with the current one
        if (colors.find(akmer) != colors.end() && ac != 255) {
            rcolor |= colors[akmer][ac];
        }
        else {
            // otherwise it is a B_{+,1} or B_{+,+} branching node => put the incoming edges to the queue and continue
            ac = bits_to_id((uint8_t) (akmer >> (kmer_bits - LOGSIGMA)).to_ulong());
            akmer = (akmer << LOGSIGMA) & mask;
            // put the actual node into the visited vector
            visited[akmer];
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                akmer ^= symbol_to_bits(base[i]);
                const uint8_t shifted_index = (uint8_t) (1 << ac);
                if (dbg_kmers.find(akmer) != dbg_kmers.end() && (dbg_kmers[akmer] & shifted_index) != 0) {
                    if ((visited[akmer] & shifted_index) == 0) {
                        kmer_queue.push_back(akmer);
                    }
                    // make sure that we won't go the same direction twice...
                    visited[akmer] |= shifted_index;
                }
                akmer ^= symbol_to_bits(base[i]);
            }
        }
    }

    return rcolor;
}


// from a given bitstring generates the appropriate kmer string
template<uint16_t KMERBITS>
string DeBrujinGraph<KMERBITS>::kmer_to_str(bitset<KMERBITS> kmer_str) {
    static const bitset<KMERBITS> mask(string(LOGSIGMA, '1'));
    stringstream ss;
    // ss << str.to_string() << endl;
    for (int i = 0; i < km; ++i) {
        uint8_t ac = (uint8_t) (kmer_str & mask).to_ulong();
        ss << bits_to_char(ac);
        kmer_str >>= LOGSIGMA;
    }
    return ss.str();
}


template<uint16_t KMERBITS>
inline bitset<MAXCOLORS> DeBrujinGraph<KMERBITS>::color_to_bitset(const Roaring &rc) {
    bitset<MAXCOLORS> acolor;
    for (uint32_t j = 0; j < C; ++j) {
        acolor[j] = rc.contains(j);
    }
    return acolor;
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::print_node(const bitset<KMERBITS> &str, uint64_t icnt, uint64_t ocnt) {
    static const bitset<KMERBITS> mask(string(kmer_bits, '1'));
    cerr << bitset<KMERBITS>(str).to_string() << endl << kmer_to_str(str) << endl;

    // print the in edges
    cout << endl << "in edges:  ";
    bitset<KMERBITS> akmer = str;
    // the last added character
    uint8_t ac = bits_to_id((uint8_t) (akmer >> (kmer_bits - LOGSIGMA)).to_ulong());
    // generate all the possible (SIGMA + 1) k-mers that could be the source
    akmer = (akmer << LOGSIGMA) & mask;
    for (int i = 0; i < SIGMA + 1; ++i) {
        akmer ^= symbol_to_bits(base[i]);
        // if (dbg_kmers.contains(akmer) && nodes_outgoing[dbg_kmers.rank(akmer)][ac]) {
        if (dbg_kmers.find(akmer) != dbg_kmers.end() && (dbg_kmers[akmer] & (1 << ac)) != 0) {
            cout << base[i] << " ";
        }
        akmer ^= symbol_to_bits(base[i]);
    }

    // print the out edges
    cout << endl << "out edges: ";
    const auto &anode = dbg_kmers[str];
    for (int i = 0; i < SIGMA + 1; ++i) {
        if (((1 << i) & anode) != 0) {
            cout << base[i] << " ";
        }
    }
    cout << endl;

    if (icnt == 1 && ocnt == 1) {
        cout << "in out, 1" << endl;
    }
    else if (icnt > 1 && ocnt <= 1) {
        cout << "branching in > 1" << endl;
    }
    else if (ocnt > 1 && icnt <= 1) {
        cout << "branching out > 1" << endl;
    }
    else if (icnt > 1 && ocnt > 1) {
        cout << "branching in and out > 1" << endl;
    }

    if (icnt == 0) {
        cout << "source" << endl;
    }
    else if (ocnt == 0) {
        cout << "sink" << endl;
    }

    cout << endl;
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
