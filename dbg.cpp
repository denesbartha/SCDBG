#include "dbg.h"

#define MAXCOLORS  95146

template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::process_read(const string &dna_str, const uint32_t color_id, bool phase_first) {
    // static T sid;
    cerr << "process " << color_id << "..." << endl;
    // read the first KMER-sized substring will be $$$..$
    // boost::multiprecision::uint256_t akmer = 0;
    bitset<KMERBITS> akmer = 0;

    if (!phase_first && !second_phase_started) {
        // nodes_outgoing.resize(dbg_kmers.size(), 0);
        // colors.resize(dbg_kmers.size());
        // for (uint64_t i = 0; i < dbg_kmers.size(); ++i) {
        //     nodes_outgoing[i] = 0;
        //
        // }
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
        // the last edge leads to $...
        dbg_kmers[akmer] |= 1;
    }
    else {
        for (uint64_t i = 0, sl = dna_str.length(); i < sl; ++i) {
            uint8_t sid = symbol_to_id(dna_str[i]);
            if (sid != 255) {
                akmer >>= LOGSIGMA;
                akmer |= shifted_sids[sid];
                if (colors.find(akmer) != colors.end()) {
                    colors[akmer][sid].add(color_id);
                }
            }
        }
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
void DeBrujinGraph<KMERBITS>::gen_succinct_dbg(const string fname) {
    cerr << "Generating Succinct De Bruijn Graph..." << endl;
    cerr << "Generating Edge list..." << endl;

    for (auto it : dbg_kmers) {
        // cout << it.first.to_string() << endl;
        for (uint8_t i = 0; i < SIGMA + 1; ++i) {
            if (((1 << i) & it.second) != 0) {
                cout << base[i];
            }
        }
    }

    // cerr << "Generating B_F list..." << endl;
    // for (auto it : dbg_kmers) {
    //     // const auto &anode = nodes_ingoing[dbg_kmers.rank(it)];
    //     int ic = indegree(it);
    //     if (ic > 0) {
    //         cout << string(ic - 1, '0') << "1";
    //     }
    //     else {
    //         cout << " ";
    //     }
    // }
    //
    // cerr << endl << "Generating B_L list..." << endl;
    // for (auto it : dbg_kmers) {
    //     const auto &anode = nodes_outgoing[dbg_kmers.rank(it)];
    //     cout << string(outdegree(anode) - 1, '0') << "1";
    // }
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::do_sampling() {
    // maximum distance without
    int max_distance = max(sampling_max_distance, (int) log2(dbg_kmers.size()));
    cerr << "Starting sampling process with max distance: " << max_distance << "..." << std::endl;
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    for (auto it = dbg_kmers.cbegin(); it != dbg_kmers.cend(); ++it) {
        uint32_t outdeg = outdegree(dbg_kmers[it->first]);
        if ((!visited[it->first] && outdeg > 1) || it == dbg_kmers.cbegin()) {
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
                        while (outdegree(dbg_kmers[akmer]) == 1) { // && indegree(akmer) == 1
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
    // uint64_t
    // uint64_t colored_branches = 0;
    // sparse_hash_map<bitset<COLORBITS>, bitset<1> > colors;
    // sparse_hash_map<bitset<KMERBITS>, uint8_t> cs;

    // for (const auto &item : dbg_kmers) {
    sparse_hash_map<bitset<MAXCOLORS>, size_t> cm;
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    sparse_hash_map<size_t, size_t> paths;
    for (auto it = dbg_kmers.cbegin(); it != dbg_kmers.cend(); ++it) {
        // bitset<kmer_bits> bb = item;
        // cout << bb.to_string() << " " << nodes_outgoing[dbg_kmers.rank(item) - 1] << endl;

        // uint64_t icnt = count_edges(item.second->in_edges);
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
                    while (outdegree(dbg_kmers[akmer]) == 1) { // && indegree(akmer) == 1
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
                }
            }
        }
        visited[it->first] = true;

        // color stats...
        // cerr << "Generate color stats..." << endl;
        if (colors.find(it->first) != colors.end()) {
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (((1 << i) & it->second) != 0) {
                    bitset<MAXCOLORS> acolor;
                    auto &rc = colors[it->first][i];
                    for (uint32_t j = 0; j < C; ++j) {
                        acolor[j] = rc.contains(j);
                    }
                    // cout << colors[it->first][i];
                    // bitset<100> bss(colors[it->first][i].toString());
                    cm[acolor]++;
                }
            }
        }
    }


    cout << "in degrees:" << endl;
    for (int i = 0; i < SIGMA + 1; ++i) {
        cout << i << ": " << in_degrees[i] << endl;
    }

    cout << "out degrees:" << endl;
    for (int i = 0; i < SIGMA + 1; ++i) {
        cout << i << ": " << out_degrees[i] << endl;
    }

    cout << fixed;
    cout << setprecision(5);
    cout << "k-mer size:\t\t\t" << (uint64_t) km << endl;
    cout << "# of nodes:\t\t\t" << num_of_nodes << endl;
    cout << "# of edges:\t\t\t" << (num_of_edges / 2) << endl;
    cout << "# of colors:\t\t\t" << (uint64_t) C << endl;
    cout << "# of color classes:\t\t" << cm.size() << endl;
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

    cout << "maximal path legth:\t\t" << paths.crbegin()->first << " " << paths.rbegin()->second << endl;

    size_t ps = 0;
    size_t os = 0;
    for (auto it = paths.cbegin(); it != paths.cend(); ++it) {
        ps += it->first * it->second;
        os += it->second;
    }
    cout << "average path legth:\t\t" << (double) ps / os << " " << paths.rbegin()->second << endl;
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
    uint8_t ac = symbol_to_id(bits_to_char((uint8_t) (pkmer >> (kmer_bits - LOGSIGMA)).to_ulong()));
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
void DeBrujinGraph<KMERBITS>::print_node(const bitset<KMERBITS> &str, uint64_t icnt, uint64_t ocnt) {
    static const bitset<KMERBITS> mask(string(kmer_bits, '1'));
    cerr << bitset<KMERBITS>(str).to_string() << endl << kmer_to_str(str) << endl;

    // print the in edges
    cout << endl << "in edges:  ";
    bitset<KMERBITS> akmer = str;
    // the last added character
    uint8_t ac = symbol_to_id(bits_to_char((uint8_t) (akmer >> (kmer_bits - LOGSIGMA)).to_ulong()));
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


// define static variables
template<uint16_t KMERBITS>
uint16_t DeBrujinGraph<KMERBITS>::kmer_bits = 0;


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
