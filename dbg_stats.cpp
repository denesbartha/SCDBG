#include "dbg_stats.h"

template<uint16_t KMERBITS>
void DeBrujinGraphStats<KMERBITS>::do_stats() {
    cerr << "Creating statistics..." << endl;

    uint64_t num_of_nodes = this->dbg_kmers.size();
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
    for (auto it = this->dbg_kmers.cbegin(); it != this->dbg_kmers.cend(); ++it) {
        uint64_t icnt = this->indegree(it->first);
        uint64_t ocnt = this->outdegree(this->dbg_kmers[it->first]);
        in_degrees[icnt]++;
        out_degrees[ocnt]++;

        // this->colors[item.second.this->colors] = 1;
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
        if ((visited.find(it->first) == visited.end() && ocnt > 1) || it == this->dbg_kmers.cbegin()) {
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (((1 << i) & it->second) != 0) {
                    size_t path_length = 1;

                    bitset<KMERBITS> akmer = it->first;
                    akmer >>= LOGSIGMA;
                    akmer |= this->shifted_sids[i];
                    while (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() &&
                           this->outdegree(this->dbg_kmers[akmer]) == 1) {
                        // && this->indegree(akmer) == 1
                        visited[akmer] = 1;
                        auto &anode = this->dbg_kmers[akmer];
                        for (uint8_t j = 0; j < SIGMA + 1; ++j) {
                            if (((1 << j) & anode) != 0) {
                                akmer >>= LOGSIGMA;
                                akmer |= this->shifted_sids[j];
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
        if (this->colors.find(it->first) != this->colors.end()) {
            if (ocnt > 1) {
                // auto current_node_color = get_color(it->first);
                Roaring prev_color;
                bool first_node = true;
                bool color_is_the_same = true;
                for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                    if (((1 << i) & it->second) != 0) {
                        auto ac = this->colors[it->first][i];
                        cm[this->color_to_bitset(ac)]++;
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
                // cm_rest[this->color_to_bitset(get_color(it->first))]++;
                for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                    if (((1 << i) & it->second) != 0) {
                        auto ac = this->colors[it->first][i];
                        cm_rest[this->color_to_bitset(ac)]++;
                    }
                }
            }
        }
        // else if (icnt > 1 || ocnt > 1) {
        //     cm_rest[this->color_to_bitset(get_color(it->first))]++;
        //
        //     if (ocnt > 1) {
        //
        //     }
        // }
    }


    cout << "in degrees:" << endl;
    for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        cout << (int) i << ": " << in_degrees[i] << endl;
    }

    cout << "out degrees:" << endl;
    for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        cout << (int) i << ": " << out_degrees[i] << endl;
    }

    cout << fixed;
    cout << setprecision(5);
    cout << "k-mer size:\t\t\t" << (uint64_t) this->km << endl;
    cout << "# of nodes:\t\t\t" << num_of_nodes << endl;
    cout << "# of edges:\t\t\t" << this->num_of_edges << endl;
    cout << "# of colors:\t\t\t" << (uint64_t) this->C << endl;
    cout << "# of color classes:\t\t" << cm.size() << endl;
    cout << "# of color classes - not stored:\t" << sparse_hash_map_difference<bitset<MAXCOLORS>, size_t>(cm_rest, cm)
         << endl;
    for (uint8_t i = 2; i <= SIGMA + 1; ++i) {
        cout << "# of B_{*," << (int) i << "} nodes, where the color is the same:\t" << branching_color_is_the_same[i]
             << "/" << branching_outdegrees[i] << endl;
    }
    // cout << "# of B_{*,+} nodes, where the outgoing color is the same as the node's color: "
    //      << color_is_the_same_as_incoming << "/" << branching_outdegrees_sum << endl;
    cout << "# of explicitly stored this->colors:\t" << this->explicitly_stored_colors << endl;
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


// from a given bitstring generates the appropriate kmer string
template<uint16_t KMERBITS>
string DeBrujinGraphStats<KMERBITS>::kmer_to_str(bitset<KMERBITS> kmer_str) {
    static const bitset<KMERBITS> mask(string(LOGSIGMA, '1'));
    stringstream ss;
    // ss << str.to_string() << endl;
    for (int i = 0; i < this->km; ++i) {
        uint8_t ac = (uint8_t) (kmer_str & mask).to_ulong();
        ss << bits_to_char(ac);
        kmer_str >>= LOGSIGMA;
    }
    return ss.str();
}


template<uint16_t KMERBITS>
void DeBrujinGraphStats<KMERBITS>::print_node(const bitset<KMERBITS> &str, uint64_t icnt, uint64_t ocnt) {
    static const bitset<KMERBITS> mask(string(this->kmer_bits, '1'));
    cerr << bitset<KMERBITS>(str).to_string() << endl << kmer_to_str(str) << endl;

    // print the in edges
    cout << endl << "in edges:  ";
    bitset<KMERBITS> akmer = str;
    // the last added character
    uint8_t ac = bits_to_id((uint8_t) (akmer >> (this->kmer_bits - LOGSIGMA)).to_ulong());
    // generate all the possible (SIGMA + 1) k-mers that could be the source
    akmer = (akmer << LOGSIGMA) & mask;
    for (int i = 0; i < SIGMA + 1; ++i) {
        akmer ^= symbol_to_bits(base[i]);
        // if (this->dbg_kmers.contains(akmer) && nodes_outgoing[this->dbg_kmers.rank(akmer)][ac]) {
        if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() && (this->dbg_kmers[akmer] & (1 << ac)) != 0) {
            cout << base[i] << " ";
        }
        akmer ^= symbol_to_bits(base[i]);
    }

    // print the out edges
    cout << endl << "out edges: ";
    const auto &anode = this->dbg_kmers[str];
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


// Returns the colour of a given node
template<uint16_t KMERBITS>
inline Roaring DeBrujinGraphStats<KMERBITS>::get_color(const bitset<KMERBITS> &pkmer) {
    static const bitset<KMERBITS> mask(string(this->kmer_bits, '1'));
    Roaring rcolor;
    deque<bitset<KMERBITS>> kmer_queue(1, pkmer);
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    while (!kmer_queue.empty()) {
        auto akmer = kmer_queue.front();
        kmer_queue.pop_front();
        std::vector<bitset<KMERBITS>> incoming_nodes;
        // step backwards until we find the color information or hit a branching node
        uint8_t ac = 255;
        while (this->colors.find(akmer) == this->colors.end() && indegree(akmer) <= 1) {
            ac = bits_to_id((uint8_t) (akmer >> (this->kmer_bits - LOGSIGMA)).to_ulong());
            // cerr << kmer_to_str(akmer) << " " << (akmer >> (kmer_bits - LOGSIGMA)).to_ulong() << endl;
            akmer = (akmer << LOGSIGMA) & mask;
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                akmer ^= symbol_to_bits(base[i]);
                if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() && (this->dbg_kmers[akmer] & (1 << ac)) != 0) {
                    // there is at most only one incoming node => break if we have found that
                    break;
                }
                akmer ^= symbol_to_bits(base[i]);
            }
        }
        // if we have found the color => merge it with the current one
        if (this->colors.find(akmer) != this->colors.end() && ac != 255) {
            rcolor |= this->colors[akmer][ac];
        }
        else {
            // otherwise it is a B_{+,1} or B_{+,+} branching node => put the incoming edges to the queue and continue
            ac = bits_to_id((uint8_t) (akmer >> (this->kmer_bits - LOGSIGMA)).to_ulong());
            akmer = (akmer << LOGSIGMA) & mask;
            // put the actual node into the visited vector
            visited[akmer];
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                akmer ^= symbol_to_bits(base[i]);
                const uint8_t shifted_index = (uint8_t) (1 << ac);
                if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end()
                    && (this->dbg_kmers[akmer] & shifted_index) != 0) {
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
