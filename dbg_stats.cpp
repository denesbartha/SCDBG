#include "dbg_stats.h"
#include <random>

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

    // sparse_hash_map<bitset<MAXCOLORS>, size_t> cm;
    // sparse_hash_map<bitset<MAXCOLORS>, size_t> cm_rest;
    // sparse_hash_map<bitset<MAXCOLORS>, bool> cm_all;
    size_t cm_color_class_sizes = 0;
    size_t cm_rest_color_class_sizes = 0;
    sparse_hash_map<size_t, uint8_t> cm;
    sparse_hash_map<size_t, uint8_t> cm_rest;
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


        // // traverse the paths
        // if ((visited.find(it->first) == visited.end() && ocnt > 1) || it == this->dbg_kmers.cbegin()) {
        //     for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        //         if (((1 << i) & it->second) != 0) {
        //             size_t path_length = 1;
        //
        //             bitset<KMERBITS> akmer = it->first;
        //             akmer >>= LOGSIGMA;
        //             akmer |= this->shifted_sids[i];
        //             while (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() &&
        //                    this->outdegree(this->dbg_kmers[akmer]) == 1) {
        //                 // && this->indegree(akmer) == 1
        //                 visited[akmer] = 1;
        //                 auto& anode = this->dbg_kmers[akmer];
        //                 for (uint8_t j = 0; j < SIGMA + 1; ++j) {
        //                     if (((1 << j) & anode) != 0) {
        //                         akmer >>= LOGSIGMA;
        //                         akmer |= this->shifted_sids[j];
        //                         ++path_length;
        //                         break;
        //                     }
        //                 }
        //             }
        //
        //             if (paths.find(path_length) != paths.end()) {
        //                 paths[path_length]++;
        //             }
        //             else {
        //                 paths[path_length] = 1;
        //             }
        //             if (path_length > maxpath_length) {
        //                 maxpath_length = path_length;
        //             }
        //         }
        //     }
        // }
        // visited[it->first] = 1;

        // color stats...
        // cerr << "Generate color stats..." << endl;
        if (this->colors.find(it->first) != this->colors.end()) {
            if (ocnt > 1) {
                // auto current_node_color = get_color(it->first);
                // Roaring prev_color;
                for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                    if (((1 << i) & it->second) != 0) {
                        auto ac = this->colors[it->first][i];
                        if (cm.find(ac) == cm.end()) {
                            cm_color_class_sizes += this->color_classes[ac].bitvector.count();
                            cm[ac];
                        }

                        // auto acolor = this->color_to_bitset(ac);
                        // cm[acolor]++;
                        // if (cm_all.find(acolor) == cm_all.end()) {
                        //     cm_color_class_sizes += ac.cardinality();
                        //     cm_all[acolor];
                        // }
                    }
                }
                branching_outdegrees[ocnt]++;
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

    for (auto it = this->dbg_kmers.cbegin(); it != this->dbg_kmers.cend(); ++it) {
        uint64_t icnt = this->indegree(it->first);
        if (this->colors.find(it->first) != this->colors.end()) {
            if (icnt > 1) {
                for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                    if (((1 << i) & it->second) != 0) {
                        auto ac = this->colors[it->first][i];
                        if (cm.find(ac) == cm.end() && cm_rest.find(ac) == cm_rest.end()) {
                            cm_rest_color_class_sizes += this->color_classes[ac].bitvector.count();
                            cm_rest[ac];
                        }
                        // auto acolor = this->color_to_bitset(ac);
                        // cm_rest[acolor]++;
                        // if (cm_all.find(acolor) == cm_all.end()) {
                        //     cm_rest_color_class_sizes += ac.cardinality();
                        //     cm_all[acolor];
                        // }
                    }
                }
            }
        }
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
    cout << "# of colors:\t\t\t" << (uint64_t) this->num_of_colors << endl;
    cout << "# of color classes:\t\t" << cm.size() << endl;
    // size_t not_stored_color_classes_cnt = sparse_hash_map_difference<bitset<MAXCOLORS>, size_t>(cm_rest, cm);
    size_t not_stored_color_classes_cnt = this->color_classes.size() - cm.size();
    cout << "# of color classes - not stored:\t" << not_stored_color_classes_cnt << endl;
    for (uint8_t i = 2; i <= SIGMA + 1; ++i) {
        cout << "# of B_{*," << (int) i << "} nodes, where the color is the same:\t" << branching_color_is_the_same[i]
             << "/" << branching_outdegrees[i] << endl;
    }
    // cout << "# of B_{*,+} nodes, where the outgoing color is the same as the node's color: "
    //      << color_is_the_same_as_incoming << "/" << branching_outdegrees_sum << endl;

    // explicitly_stored_colors - shows that how many places did we store labels
    cout << "# of explicitly stored color labels:\t" << this->explicitly_stored_colors << endl;
    cout << "color class sizes - stored:\t" << cm_color_class_sizes << "\t average:"
         << (cm_color_class_sizes / (float) cm.size()) << endl;
    cout << "color class sizes - not stored:\t" << cm_rest_color_class_sizes << "\t average:"
         << (cm_rest_color_class_sizes / (float) not_stored_color_classes_cnt) << endl;
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
    cout << "average path legth:\t\t" << (double) ps / os << endl << endl;


    static const bitset<KMERBITS> mask(string(this->kmer_bits, '1'));
    {
        sparse_hash_map<size_t, size_t> trees;
        sparse_hash_map<size_t, size_t> leaves_avg;
        uint64_t maxtree_size = 0;
        uint64_t maxleaves_cnt = 0;

        paths.clear();
        maxpath_length = 0;

        if (this->dbg_kmers_sorted.size() == 0) {
            cerr << "sorting..." << endl;
            this->sort_dbg();

            // for (size_t i = 0; i < this->dbg_kmers_sorted.size(); ++i) {
            //     auto& anode = this->dbg_kmers_sorted[i].second;
            //     for (uint8_t j = 0; j < SIGMA + 1; ++j) {
            //         if (((1 << j) & anode) != 0) {
            //             cout << i << " " << kmer_to_str(this->dbg_kmers_sorted[i].first) << " " << base[j] << endl;
            //         }
            //     }
            // }
            // cout << endl;
        }


        cerr << "Gen runtime stats..." << endl;
        // Seed with a real random value, if available
        // random_device rd;
        // mt19937 engine(rd());
        // uniform_int_distribution<size_t> uniform_dist(0, this->dbg_kmers_sorted.size() - 1);
        cerr << "gen stats on random edges..." << endl;
        for (int i = 0; i < 100000; ++i) {

            bitset<KMERBITS> kmer;
            size_t index = rand() % this->dbg_kmers_sorted.size(); // uniform_dist(engine); //

            kmer = this->dbg_kmers_sorted[index].first;
            uint8_t indegree = this->indegree(kmer);
            if (indegree == 0) continue;

            uniform_int_distribution<> uniform_dist2(0, indegree - 1);
            uint8_t j = rand() % indegree; // uniform_dist2(engine);


            // cout << "first :" << this->kmer_to_str(kmer) << endl;

            uint8_t ac = bits_to_id((uint8_t) (kmer >> (this->kmer_bits - LOGSIGMA)).to_ulong());
            kmer = (kmer << LOGSIGMA) & mask;
            for (uint8_t jj = 0; jj < SIGMA + 1; ++jj) {
                kmer ^= symbol_to_bits(base[jj]);
                // cout << this->kmer_to_str(kmer) << endl;
                if (this->dbg_kmers.find(kmer) != this->dbg_kmers.end() && (this->dbg_kmers[kmer] & (1 << ac)) != 0) {
                    if (j-- == 0) {

                        sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
                        auto p2 = traverse(kmer, false, visited);
                        size_t tree_size = get<0>(p2);
                        size_t leaf_cnt = get<1>(p2);
                        size_t path_size = get<2>(p2);
                        // cout << this->kmer_to_str(kmer) << " " << tree_size << " " << leaf_cnt << " " << path_size << endl;
                        if (trees.find(tree_size) != trees.end()) {
                            trees[tree_size]++;
                        }
                        else {
                            trees[tree_size] = 1;
                        }
                        if (tree_size > maxtree_size) {
                            maxtree_size = tree_size;
                        }

                        if (leaves_avg.find(leaf_cnt) != leaves_avg.end()) {
                            leaves_avg[leaf_cnt]++;
                        }
                        else {
                            leaves_avg[leaf_cnt] = 1;
                        }
                        if (leaf_cnt > maxleaves_cnt) {
                            maxleaves_cnt = leaf_cnt;
                        }

                        if (paths.find(path_size) != paths.end()) {
                            paths[path_size]++;
                        }
                        else {
                            paths[path_size] = 1;
                        }
                        if (path_size > maxpath_length) {
                            maxpath_length = path_size;
                        }
                        break;
                    }
                }
                kmer ^= symbol_to_bits(base[jj]);
            }


            if (i % 1000 == 0) {
                cerr << i << " ";

                cout << "maximal tree size:\t\t" << maxtree_size << " " << trees[maxtree_size] << endl;
                size_t ps = 0;
                size_t os = 0;
                for (auto it = trees.cbegin(); it != trees.cend(); ++it) {
                    ps += it->first * it->second;
                    os += it->second;
                }
                cout << "average tree size:\t\t" << (double) ps / os << endl;


                cout << "maximal leaf cnt:\t\t" << maxleaves_cnt << " " << leaves_avg[maxleaves_cnt] << endl;
                ps = 0;
                os = 0;
                for (auto it = leaves_avg.cbegin(); it != leaves_avg.cend(); ++it) {
                    ps += it->first * it->second;
                    os += it->second;
                }
                cout << "average leaf cnt:\t\t" << (double) ps / os << endl;


                cout << "maximal path size:\t\t" << maxpath_length << " " << paths[maxpath_length] << endl;
                ps = 0;
                os = 0;
                for (auto it = paths.cbegin(); it != paths.cend(); ++it) {
                    ps += it->first * it->second;
                    os += it->second;
                }
                cout << "average path size:\t\t" << (double) ps / os << endl;
            }
        }
        cerr << endl;



    }

    // cout << endl;
    // this->sort_dbg();
    // for (auto it = this->dbg_kmers_sorted.cbegin(); it != this->dbg_kmers_sorted.cend(); ++it) {
    //     cout << kmer_to_str(it->first) << endl;
    // }
}


template<uint16_t KMERBITS>
tuple<size_t, size_t, size_t> DeBrujinGraphStats<KMERBITS>::traverse(bitset<KMERBITS> pkmer, bool path,
                                                            sparse_hash_map<bitset<KMERBITS>, uint8_t>& visited) {
    static const bitset<KMERBITS> mask(string(this->kmer_bits, '1'));
    size_t size = 1;

    bitset<KMERBITS> akmer = pkmer;
    size_t leaves = 0;
    size_t path_cnt = 0;

    visited[akmer];
    // cout << this->kmer_to_str(akmer) << " " << (int) this->indegree(akmer) << " " << (this->colors.find(akmer) == this->colors.end()) << endl;
    while (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() && this->indegree(akmer) == 1
           && this->colors.find(akmer) == this->colors.end()) { //&& this->outdegree(akmer) == 1) {

        uint8_t ac = bits_to_id((uint8_t) (akmer >> (this->kmer_bits - LOGSIGMA)).to_ulong());
        akmer = (akmer << LOGSIGMA) & mask;
        for (uint8_t j = 0; j < SIGMA + 1; ++j) {
            akmer ^= symbol_to_bits(base[j]);
            if (visited.find(akmer) != visited.end()) {
                break;
            }
            if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() && (this->dbg_kmers[akmer] & (1 << ac)) != 0) {
                ++size;
                // cout << this->kmer_to_str(akmer) << " " << (int) this->indegree(akmer) << " " << (this->colors.find(akmer) == this->colors.end()) << endl;
                break;
            }
            akmer ^= symbol_to_bits(base[j]);
        }

        if (visited.find(akmer) != visited.end()) {
            break;
        }

        visited[akmer];
    }
    path_cnt = size;
    if (path) {
        return tuple<size_t, size_t, size_t>(size, 1, path_cnt);
    }
    else {
        visited[akmer];
        if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() && this->indegree(akmer) > 1 &&
            this->colors.find(akmer) == this->colors.end()) { // visited.find(akmer) == visited.end()

            uint8_t ac = bits_to_id((uint8_t) (akmer >> (this->kmer_bits - LOGSIGMA)).to_ulong());
            akmer = (akmer << LOGSIGMA) & mask;
            for (uint8_t j = 0; j < SIGMA + 1; ++j) {
                akmer ^= symbol_to_bits(base[j]);
                if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() && (this->dbg_kmers[akmer] & (1 << ac)) != 0) {
                    if (visited.find(akmer) == visited.end()) {
                        auto p = traverse(akmer, false, visited);
                        size += get<0>(p);
                        leaves += get<1>(p);
                    }
                }
                akmer ^= symbol_to_bits(base[j]);
            }
        }
        else if (this->colors.find(akmer) != this->colors.end()) {
            leaves++;
        }

        return tuple<size_t, size_t, size_t>(size, leaves, path_cnt);
    }
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
void DeBrujinGraphStats<KMERBITS>::print_node(const bitset<KMERBITS>& str, uint64_t icnt, uint64_t ocnt) {
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
    const auto& anode = this->dbg_kmers[str];
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
inline void DeBrujinGraphStats<KMERBITS>::get_color(const bitset<KMERBITS>& pkmer) {
    static const bitset<KMERBITS> mask(string(this->kmer_bits, '1'));
    // Roaring rcolor;
    size_t leaf_cnt = 0;
    size_t edge_cnt = 0;
    deque<bitset<KMERBITS>> kmer_queue(1, pkmer);
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    while (!kmer_queue.empty()) {
        auto akmer = kmer_queue.front();
        kmer_queue.pop_front();
        std::vector<bitset<KMERBITS>> incoming_nodes;
        // step backwards until we find the color information or hit a branching node
        visited[akmer];
        uint8_t ac = 255;
        while (this->colors.find(akmer) == this->colors.end() && indegree(akmer) <= 1) {
            visited[akmer];
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
            // rcolor |= this->colors[akmer][ac];
            leaf_cnt++;
        }
        else {
            // otherwise it is a B_{+,1} or B_{+,+} branching node => put the incoming edges to the queue and continue
            ac = bits_to_id((uint8_t) (akmer >> (this->kmer_bits - LOGSIGMA)).to_ulong());
            akmer = (akmer << LOGSIGMA) & mask;
            // put the actual node into the visited vector
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                akmer ^= symbol_to_bits(base[i]);
                const uint8_t shifted_index = (uint8_t) (1 << ac);
                if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() && (this->dbg_kmers[akmer] & shifted_index) != 0) {
                    if (visited.find(akmer) == visited.end()) { // & shifted_index) == 0
                        kmer_queue.push_back(akmer);
                        edge_cnt++;
                    }
                    // make sure that we won't go the same direction twice...
                    // visited[akmer] |= shifted_index;
                }
                akmer ^= symbol_to_bits(base[i]);
            }
        }
    }

    // return rcolor;
}
