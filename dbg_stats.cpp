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
    uint64_t maxtree_size = 0;
    uint64_t maxleaves_cnt = 0;
    uint64_t maxpath_length = 0;
    // uint64_t branching_color_is_the_same[SIGMA + 2] = {};
    // uint64_t branching_outdegrees[SIGMA + 2] = {};
    // uint64_t color_is_the_same_as_incoming = 0;
    // uint64_t branching_outdegrees_sum = 0;

    // sparse_hash_map<bitset<MAXCOLORS>, size_t> cm;
    // sparse_hash_map<bitset<MAXCOLORS>, size_t> cm_rest;
    // sparse_hash_map<bitset<MAXCOLORS>, bool> cm_all;
    size_t cm_color_class_sizes = 0;
    size_t cm_rest_color_class_sizes = 0;
    sparse_hash_map<size_t, uint8_t> cm;
    sparse_hash_map<size_t, uint8_t> cm_rest;
    // sparse_hash_map<bitset<KMERBITS>, size_t> visited;

    // traverse the trees
    std::function<pair<size_t, size_t>(bitset<KMERBITS>, uint8_t, bool)> traverse;
    traverse = [this, &traverse](bitset<KMERBITS> pkmer, uint8_t i, bool path) -> pair<size_t, size_t> {
        size_t size = 1;

        bitset<KMERBITS> akmer = pkmer;
        size_t leaves = 0;
        akmer >>= LOGSIGMA;
        akmer |= this->shifted_sids[i];
        while (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() &&
               this->outdegree(this->dbg_kmers[akmer]) == 1 && this->indegree(akmer) == 1) {
            auto& anode = this->dbg_kmers[akmer];
            for (uint8_t j = 0; j < SIGMA + 1; ++j) {
                if (((1 << j) & anode) != 0) {
                    akmer >>= LOGSIGMA;
                    akmer |= this->shifted_sids[j];
                    ++size;
                    break;
                }
            }
        }
        if (path) {
            return pair<size_t, size_t>(size, 1);
        }
        else {
            if (this->dbg_kmers.find(akmer) != this->dbg_kmers.end() &&
                this->indegree(akmer) == 1 && this->outdegree(this->dbg_kmers[akmer]) > 1) {

                for (uint8_t j = 0; j < SIGMA + 1; ++j) {
                    if (((1 << j) & this->dbg_kmers[akmer]) != 0) {
                        auto p = traverse(akmer, j, false);
                        size += p.first;
                        leaves += p.second;
                    }
                }
            }
            if (this->dbg_kmers.find(akmer) == this->dbg_kmers.end() || this->indegree(akmer) > 1) {
                leaves++;
            }
            return pair<size_t, size_t>(size, leaves);
        }
    };

    // sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    sparse_hash_map<size_t, size_t> trees;
    sparse_hash_map<size_t, size_t> leaves_avg;
    sparse_hash_map<size_t, size_t> paths_avg;
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


        if (icnt > 1 || it == this->dbg_kmers.cbegin()) {
            uint8_t edges = this->dbg_kmers[it->first];
            for (uint8_t j = 0; j < SIGMA + 1; ++j) {
                if (((1 << j) & edges) != 0) {
                    // traverse like a tree...
                    auto p = traverse(it->first, j, false);
                    size_t tree_size = p.first;
                    size_t leaf_cnt = p.second;
                    // cout << this->kmer_to_str(it->first) << base[j] << " " << tree_size << " " << leaf_cnt << endl;
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
                }
            }
        }

        if (icnt > 1 || ocnt > 1 || it == this->dbg_kmers.cbegin()) {
            uint8_t edges = this->dbg_kmers[it->first];
            for (uint8_t j = 0; j < SIGMA + 1; ++j) {
                if (((1 << j) & edges) != 0) {
                    // traverse the paths only...
                    size_t path_length = traverse(it->first, j, true).first;
                    if (paths_avg.find(path_length) != paths_avg.end()) {
                        paths_avg[path_length]++;
                    }
                    else {
                        paths_avg[path_length] = 1;
                    }
                    if (path_length > maxpath_length) {
                        maxpath_length = path_length;
                    }

                    // cout << this->kmer_to_str(it->first) << base[j] << " " << path_length << endl;
                }
            }
        }

        // color stats...
        // cerr << "Generate color stats..." << endl;
        if (this->colors.find(it->first) != this->colors.end()) {
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                auto ac = this->colors[it->first][i];
                if (((1 << i) & it->second) != 0 && ac != 0) {

                    bitset<KMERBITS> akmer = it->first;
                    akmer >>= LOGSIGMA;
                    akmer |= this->shifted_sids[i];

                    if (i == 0 || this->indegree(akmer) > 1) {
                        if (cm.find(ac) == cm.end()) {
                            cm_color_class_sizes += this->color_classes[ac].bitvector.count();
                            cm[ac];
                        }

                        // branching_outdegrees[ocnt]++;
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

    for (auto it = this->dbg_kmers.cbegin(); it != this->dbg_kmers.cend(); ++it) {
        if (this->colors.find(it->first) != this->colors.end()) {
            for (uint8_t i = 1; i < SIGMA + 1; ++i) {
                auto ac = this->colors[it->first][i];
                if (((1 << i) & it->second) != 0 && ac != 0) {

                    bitset<KMERBITS> akmer = it->first;
                    akmer >>= LOGSIGMA;
                    akmer |= this->shifted_sids[i];

                    if (this->indegree(akmer) == 1) {
                        if (cm.find(ac) == cm.end() && cm_rest.find(ac) == cm_rest.end()) {
                            cm_rest_color_class_sizes += this->color_classes[ac].bitvector.count();
                            cm_rest[ac];
                            // cout << this->color_classes[ac].bitvector.to_string() << endl;
                            // cout << this->kmer_to_str(it->first) << endl;
                        }

                        // branching_outdegrees[ocnt]++;
                    }
                }
            }
        }
        // // uint64_t icnt = this->indegree(it->first);
        // uint64_t ocnt = this->outdegree(this->dbg_kmers[it->first]);
        // if (this->colors.find(it->first) != this->colors.end()) {
        //     if (ocnt > 1) {
        //         for (uint8_t i = 0; i < SIGMA + 1; ++i) {
        //             if (((1 << i) & it->second) != 0) {
        //                 auto ac = this->colors[it->first][i];
        //                 if (cm.find(ac) == cm.end() && cm_rest.find(ac) == cm_rest.end()) {
        //                     cm_rest_color_class_sizes += this->color_classes[ac].bitvector.count();
        //                     cm_rest[ac];
        //                 }
        //                 // auto acolor = this->color_to_bitset(ac);
        //                 // cm_rest[acolor]++;
        //                 // if (cm_all.find(acolor) == cm_all.end()) {
        //                 //     cm_rest_color_class_sizes += ac.cardinality();
        //                 //     cm_all[acolor];
        //                 // }
        //             }
        //         }
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
    cout << "# of colors:\t\t\t" << (uint64_t) this->num_of_colors << endl;
    cout << "# of color classes:\t\t" << cm.size() << " rest size: " << cm_rest.size() << endl;
    // size_t not_stored_color_classes_cnt = sparse_hash_map_difference<bitset<MAXCOLORS>, size_t>(cm_rest, cm);
    size_t not_stored_color_classes_cnt = this->color_classes.size() - cm.size();
    cout << "# of color classes - not stored:\t" << not_stored_color_classes_cnt << endl;
    // for (uint8_t i = 2; i <= SIGMA + 1; ++i) {
    //     cout << "# of B_{*," << (int) i << "} nodes, where the color is the same:\t" << branching_color_is_the_same[i]
    //          << "/" << branching_outdegrees[i] << endl;
    // }
    // cout << "# of B_{*,+} nodes, where the outgoing color is the same as the node's color: "
    //      << color_is_the_same_as_incoming << "/" << branching_outdegrees_sum << endl;

    // explicitly_stored_labels - shows that how many places did we store labels
    cout << "# of explicitly stored color labels:\t" << this->explicitly_stored_labels << endl;
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


    cout << "maximal tree size:\t\t" << maxtree_size << " " << trees[maxtree_size] << endl;
    size_t ps = 0;
    size_t os = 0;
    for (auto it = trees.cbegin(); it != trees.cend(); ++it) {
        ps += it->first * it->second;
        os += it->second;
    }
    cout << "average tree size:\t\t" << (double) ps / os << endl;


    cout << "maximal # of leaf nodes:\t\t" << maxleaves_cnt << " " << leaves_avg[maxleaves_cnt] << endl;
    ps = 0;
    os = 0;
    for (auto it = leaves_avg.cbegin(); it != leaves_avg.cend(); ++it) {
        ps += it->first * it->second;
        os += it->second;
    }
    cout << "average # of leaves:\t\t" << (double) ps / os << endl;


    cout << "maximal path length:\t\t" << maxpath_length << " " << paths_avg[maxpath_length] << endl;
    ps = 0;
    os = 0;
    for (auto it = paths_avg.cbegin(); it != paths_avg.cend(); ++it) {
        ps += it->first * it->second;
        os += it->second;
    }
    cout << "average path length:\t\t" << (double) ps / os << endl;


    // {
    //     trees.clear();
    //     leaves_avg.clear();
    //     uint64_t maxtree_size = 0;
    //     uint64_t maxleaves_cnt = 0;
    //
    //     cerr << "sorting..." << endl;
    //     this->sort_dbg();
    //     cerr << "gen stats on random edges..." << endl;
    //     for (int i = 0; i < 1000; ++i) {
    //         auto p = this->dbg_kmers_sorted[rand() % this->dbg_kmers_sorted.size()];
    //         auto kmer = p.first;
    //         uint8_t edge = p.second;
    //         uint8_t outdeg = this->outdegree(edge);
    //         int j = rand() % outdeg;
    //
    //
    //         for (uint8_t jj = 0; jj < SIGMA + 1; ++jj) {
    //             if (((1 << jj) & edge) != 0) {
    //                 if (j-- > 0) continue;
    //
    //                 auto p2 = traverse(kmer, jj);
    //                 size_t path_length = p2.first;
    //                 size_t branching_cnt = p2.second;
    //                 // cout << this->kmer_to_str(it->first) << base[j] << " " << path_length << " " << branching_cnt << endl;
    //                 if (trees.find(path_length) != trees.end()) {
    //                     trees[path_length]++;
    //                 }
    //                 else {
    //                     trees[path_length] = 1;
    //                 }
    //                 if (path_length > maxtree_size) {
    //                     maxtree_size = path_length;
    //                 }
    //
    //                 if (leaves_avg.find(branching_cnt) != leaves_avg.end()) {
    //                     leaves_avg[branching_cnt]++;
    //                 }
    //                 else {
    //                     leaves_avg[branching_cnt] = 1;
    //                 }
    //                 if (branching_cnt > maxleaves_cnt) {
    //                     maxleaves_cnt = branching_cnt;
    //                 }
    //                 break;
    //             }
    //         }
    //
    //         if (i % 1000 == 0) {
    //             cerr << i << " ";
    //         }
    //     }
    //     cerr << endl;
    //
    //
    //     cout << "maximal tree size:\t\t" << maxtree_size << " " << trees[maxtree_size] << endl;
    //     size_t ps = 0;
    //     size_t os = 0;
    //     for (auto it = trees.cbegin(); it != trees.cend(); ++it) {
    //         ps += it->first * it->second;
    //         os += it->second;
    //     }
    //     cout << "average tree size:\t\t" << (double) ps / os << endl;
    //
    //
    //     cout << "maximal branching node cnt:\t\t" << maxleaves_cnt << " " << leaves_avg[maxleaves_cnt] << endl;
    //     ps = 0;
    //     os = 0;
    //     for (auto it = leaves_avg.cbegin(); it != leaves_avg.cend(); ++it) {
    //         ps += it->first * it->second;
    //         os += it->second;
    //     }
    //     cout << "average branching cnt:\t\t" << (double) ps / os << endl;
    // }

    // cout << endl;
    // this->sort_dbg();
    // for (auto it = this->dbg_kmers_sorted.cbegin(); it != this->dbg_kmers_sorted.cend(); ++it) {
    //     cout << kmer_to_str(it->first) << endl;
    // }
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
inline Roaring DeBrujinGraphStats<KMERBITS>::get_color(const bitset<KMERBITS>& pkmer) {
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
