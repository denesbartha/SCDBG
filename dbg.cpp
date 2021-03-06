#include <deque>
#include <sstream>
#include <algorithm>
#include <sdsl/bit_vectors.hpp>
#include <functional>
#include <stxxl.h>

#include "dbg.h"


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::process_read(const string& dna_str, const uint32_t color_id, bool phase_first) {
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
                    // colors[pkmer][sid].add(color_id);
                    add_color(colors[pkmer][sid], color_id);
                }
                pkmer = akmer;
            }
        }
        if (colors.find(pkmer) != colors.end()) {
            // colors[pkmer][0].add(color_id);
            add_color(colors[pkmer][0], color_id);
        }
    }
}


template<uint16_t KMERBITS>
inline void DeBrujinGraph<KMERBITS>::add_color(size_t& kmer_color_hash, const uint32_t color_id) {
    // if this is a new value
    if (kmer_color_hash == 0) {
        bitset<MAXCOLORS> bitvector;
        bitvector.set(color_id);
        kmer_color_hash = add_color_class(bitvector);
    }
    else if (!color_classes[kmer_color_hash].bitvector[color_id]) {
        // if the color class is no longer exists => delete it
        if (--color_classes[kmer_color_hash].cnt == 0) {
            bitset<MAXCOLORS> bitvector = color_classes[kmer_color_hash].bitvector;
            color_classes.erase(kmer_color_hash);
            set_bits -= bitvector.count();

            bitvector.set(color_id);
            kmer_color_hash = add_color_class(bitvector);
        }
        else {
            bitset<MAXCOLORS> bitvector = color_classes[kmer_color_hash].bitvector;
            bitvector.set(color_id);
            kmer_color_hash = add_color_class(bitvector);
        }
    }
}


template<uint16_t KMERBITS>
inline size_t DeBrujinGraph<KMERBITS>::add_color_class(const bitset<MAXCOLORS>& bitvector) {
    size_t hashv = hash_color_class(bitvector);
    while (color_classes.find(hashv) != color_classes.end() && color_classes[hashv].bitvector != bitvector) {
        hashv = hash_int(hashv + 1);
    }
    if (color_classes.find(hashv) == color_classes.end()) {
        color_classes[hashv].bitvector = bitvector; // insert(hashv, color_class_t(bitvector)); //
        set_bits += bitvector.count();
    }
    color_classes[hashv].cnt++;
    return hashv;
}


template<uint16_t KMERBITS>
inline void DeBrujinGraph<KMERBITS>::add_new_node(const bitset<KMERBITS>& akmer, bool new_node, uint8_t pc) {
    static bitset<KMERBITS> pkmer;
    if (!new_node) {
        dbg_kmers[pkmer] |= 1 << pc;
    }

    // Note that if there is no akmer in dbg_kmers, it will get allocated
    dbg_kmers[akmer];

    pkmer = akmer;
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::sort_dbg() {
    cerr << "Sorting edge list..." << endl;

    dbg_kmers_sorted.resize(dbg_kmers.size());
    size_t i = 0;
    for (auto it = dbg_kmers.cbegin(); it != dbg_kmers.end(); ++it, ++i) {
        dbg_kmers_sorted[i] = make_pair(it->first, it->second);

    }

    // the storting process uses 64GB ram
    // TODO: add this as an argument of the program
    stxxl::sort(dbg_kmers_sorted.begin(), dbg_kmers_sorted.end(), compare_bit_vector(kmer_bits), 64 * 1024UL * 1024UL * 1024UL);

    // std::sort(dbg_kmers_sorted.begin(), dbg_kmers_sorted.end(),
    //           [this](const pair<bitset<KMERBITS>, uint8_t>& a, const pair<bitset<KMERBITS>, uint8_t>& b) -> bool {
    //               for (int i = kmer_bits - 1; i >= 0; --i) {
    //                   if (a.first[i] ^ b.first[i]) {
    //                       return b.first[i];
    //                   }
    //               }
    //               return false;
    //           });
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::save_edge_list(ofstream& f) {
    cerr << "Saving edge list to file..." << endl;

    char data_buffer[LOGSIGMA * 8000] = {};
    size_t buffer_index = 0;
    for (size_t i = 0; i < dbg_kmers_sorted.size(); ++i) {
        for (uint8_t j = 0; j < SIGMA + 1; ++j) {
            if (((1 << j) & dbg_kmers_sorted[i].second) != 0) {
                // cout << base[j];
                uint8_t ibits = id_to_bits(j);
                for (uint8_t k = 0; k < LOGSIGMA; ++k) {
                    if ((1 << k) & ibits) {
                        data_buffer[buffer_index / 8] |= (1 << (buffer_index % 8));
                    }
                    if ((++buffer_index / 8) >= sizeof(data_buffer)) {
                        // save_data(f, sizeof(data_buffer));
                        buffer_to_file(f, data_buffer, sizeof(data_buffer), buffer_index);
                    }
                }
            }
        }
    }
    if (buffer_index > 0) {
        // f.write(data_buffer, buffer_index);
        buffer_to_file(f, data_buffer, divide_and_to_upper(buffer_index, 8));
    }
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::save_table_F(ofstream& f) {
    static const bitset<KMERBITS> mask(string(LOGSIGMA, '1'));
    cerr << "Saving table F to file..." << endl;

    array<size_t, SIGMA> F = {};
    std::fill(F.begin(), F.end(), num_of_edges);
    uint8_t sigma_id = 0;
    uint8_t pc = 0;
    // generate F table - start from index 1, skip the $, since there is only one $$$..$ node
    for (size_t i = 1; i < dbg_kmers_sorted.size() && sigma_id < SIGMA; ++i) {
        uint8_t ac = bits_to_id((uint8_t) (dbg_kmers_sorted[i].first >> (kmer_bits - LOGSIGMA)).to_ulong());
        if (ac != pc) {
            F[sigma_id] = i;
            sigma_id++;
            pc = ac;
        }
    }

    // write table F to the file
    for (uint8_t i = 0; i < SIGMA; ++i) {
        f.write(reinterpret_cast<const char *>(&F[i]), sizeof(F[i]));
    }
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::gen_succinct_dbg(const string& fname) {
    cerr << "Generating Succinct De Bruijn Graph..." << endl;

    sort_dbg();

    ofstream f(fname + ".dbg", ios::out | ios::binary);
    // write the size of the DBG to the file
    f.write(reinterpret_cast<const char *>(&num_of_edges), sizeof(num_of_edges));

    save_edge_list(f);
    save_table_F(f);

    // lambda function for saving B_L and B_F
    auto save_bin_list = [&f, this](auto fn) {
        char data_buffer[LOGSIGMA * 8000] = {};
        size_t buffer_index = 0;
        for (auto it = dbg_kmers_sorted.cbegin(); it != dbg_kmers_sorted.cend(); ++it) {
            uint8_t ic = fn(it->first);
            if (ic > 0) {
                buffer_index += ic - 1;
                if (buffer_index >= sizeof(data_buffer)) {
                    buffer_to_file(f, data_buffer, sizeof(data_buffer), buffer_index, false);
                    buffer_index = buffer_index % sizeof(data_buffer);
                }
                data_buffer[buffer_index / 8] |= (1 << (buffer_index % 8));
                ++buffer_index;
            }
        }
        if (buffer_index > 0) {
            buffer_to_file(f, data_buffer, divide_and_to_upper(buffer_index, 8));
        }
    };

    cerr << "Saving B_L list..." << endl;
    save_bin_list([this](auto kmer) { return outdegree(dbg_kmers[kmer]); });
    cerr << endl;

    cerr << "Saving B_F list..." << endl;
    save_bin_list([this](auto kmer) { return indegree(kmer); });

    f.close();

    save_colors(fname);
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::save_colors(const string& fname) {
    auto flip = [](const std::pair<uint64_t, color_class_t>& p) {
        return std::pair<size_t, bitset<MAXCOLORS>>(p.second.cnt, p.second.bitvector);
    };

    cerr << "Sorting color classes..." << endl;
    multimap<size_t, bitset<MAXCOLORS>> ordered_cm;  // flip_map(color_classes);
    std::transform(color_classes.begin(), color_classes.end(), std::inserter(ordered_cm, ordered_cm.begin()), flip);
    // color_classes.clear();

    sparse_hash_map<bitset<MAXCOLORS>, pair<size_t, uint8_t>> color_class_order;
    save_color_classes(fname, ordered_cm, color_class_order, set_bits);
    ordered_cm.clear();

    // create file for storing the bit vectors
    ofstream f(fname + ".sbv", ios::out | ios::binary);
    save_store_vector(f);
    f.close();

    cerr << "Saving label bit vector..." << endl;
    save_color_bit_vector(fname + ".lbv", color_class_order, false);

    cerr << "Saving boundary bit vector..." << endl;
    save_color_bit_vector(fname + ".bbv", color_class_order, true);
}


template<uint16_t KMERBITS>
void
DeBrujinGraph<KMERBITS>::save_color_classes(const string& fname, const multimap<size_t, bitset<MAXCOLORS>>& ordered_cm,
                                            sparse_hash_map<bitset<MAXCOLORS>, pair<size_t, uint8_t>>& color_class_order,
                                            size_t set_bits) {
    cerr << "Saving color classes..." << endl;
    auto vector_builder = new sdsl::sd_vector_builder(ordered_cm.size() * num_of_colors, set_bits);
    size_t i = 0, index, border = 2;
    uint8_t index_length = 1;
    for (auto it = ordered_cm.rbegin(); it != ordered_cm.rend(); ++it, ++i) {
        for (size_t j = 0; j < num_of_colors; ++j) {
            if (it->second[j]) {
                vector_builder->set(i * num_of_colors + j);
            }
        }
        color_class_order[it->second] = pair<size_t, uint8_t>(index, index_length);
        if (++index >= border) {
            index = 0;
            index_length++;
            border *= 2;
        }
    }

    sdsl::sd_vector<> color_matrix(*vector_builder);
    sdsl::store_to_file(color_matrix, fname + ".color_classes");
    delete vector_builder;
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::save_store_vector(ostream& f) {
    bool sampling = sampling_max_distance > 0;
    // write to file true iff there was sampling, false otherwise
    f.write(reinterpret_cast<const char *>(&sampling), sizeof(sampling));
    // if there was no sampling => we don't need to store extra information => skip this saving part
    if (!sampling) {
        cerr << "Skip saving store vector..." << endl;
        return;
    }

    cerr << "Saving store vector..." << endl;
    char data_buffer[LOGSIGMA * 8000] = {};
    size_t buffer_index = 0;
    for (auto it = dbg_kmers_sorted.cbegin(); it != dbg_kmers_sorted.cend(); ++it) {
        // if we are on a branching node (B_{1,+} or B_{+,+})
        if (colors.find(it->first) != colors.end() && outdegree(dbg_kmers[it->first]) > 1) {
            // go through the outgoing edges
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (((1 << i) & it->second) != 0) {
                    data_buffer[buffer_index / 8] |= 1 << (buffer_index % 8);
                    if ((++buffer_index / 8) >= sizeof(data_buffer)) {
                        buffer_to_file(f, data_buffer, sizeof(data_buffer), buffer_index);
                    }
                }
            }
        }
        else {
            if ((++buffer_index / 8) >= sizeof(data_buffer)) {
                buffer_to_file(f, data_buffer, sizeof(data_buffer), buffer_index);
            }
        }
    }
    if (buffer_index > 0) {
        buffer_to_file(f, data_buffer, divide_and_to_upper(buffer_index, 8));
    }
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::save_color_bit_vector(const string& fname,
                                                    sparse_hash_map<bitset<MAXCOLORS>, pair<size_t, uint8_t>>&
                                                    color_class_order,
                                                    bool boundary) {
    size_t total_bits = 0, set_bits = 0;
    for (auto it = dbg_kmers_sorted.cbegin(); it != dbg_kmers_sorted.cend(); ++it) {
        // if we are on a branching node (B_{1,+} or B_{+,+})
        if (colors.find(it->first) != colors.end()) {
            // go through the outgoing edges
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (((1 << i) & it->second) != 0) {
                    auto ac = colors[it->first][i];
                    // if the color was really used...
                    if (ac != 0) {
                        auto order_data = color_class_order[color_classes[ac].bitvector];
                        auto number_of_bits = order_data.second;
                        total_bits += number_of_bits;
                        if (boundary) {
                            set_bits++;
                        }
                        else {
                            size_t data = order_data.first;
                            for (uint8_t j = 0; j < number_of_bits; ++j) {
                                if (((1 << j) & data) != 0) {
                                    set_bits++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    auto vector_builder = new sdsl::sd_vector_builder(total_bits, set_bits);
    size_t index = 0;
    for (auto it = dbg_kmers_sorted.cbegin(); it != dbg_kmers_sorted.cend(); ++it) {
        // if we are on a branching node (B_{1,+} or B_{+,+})
        if (colors.find(it->first) != colors.end()) {
            // go through the outgoing edges
            for (uint8_t i = 0; i < SIGMA + 1; ++i) {
                if (((1 << i) & it->second) != 0) {
                    auto ac = colors[it->first][i];
                    // if the color was really used...
                    if (ac != 0) {
                        auto order_data = color_class_order[color_classes[ac].bitvector];
                        auto number_of_bits = order_data.second;
                        if (boundary) {
                            vector_builder->set(index + number_of_bits - 1);
                        }
                        else {
                            size_t data = order_data.first;
                            for (uint8_t j = 0; j < number_of_bits; ++j) {
                                if (((1 << j) & data) != 0) {
                                    vector_builder->set(index + j);
                                }
                            }
                        }
                        index += number_of_bits;
                    }
                }
            }
        }
    }

    sdsl::sd_vector<> color_bit_vector(*vector_builder);
    sdsl::store_to_file(color_bit_vector, fname);
    delete vector_builder;
}


template<uint16_t KMERBITS>
void DeBrujinGraph<KMERBITS>::do_sampling() {
    // maximum distance without
    uint32_t max_distance = max(sampling_max_distance, (uint32_t) log2(dbg_kmers.size()));
    cerr << "Starting sampling process with max distance: " << max_distance << "..." << std::endl;
    sparse_hash_map<bitset<KMERBITS>, uint8_t> visited;
    for (auto it = dbg_kmers.cbegin(); it != dbg_kmers.cend(); ++it) {
        uint32_t outdeg = outdegree(dbg_kmers[it->first]);
        num_of_edges += outdeg;
        if (!visited[it->first] && ((outdeg > 1))) { // || indegree(it->first) > 1   || it == dbg_kmers.cbegin()
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
                            auto& anode = dbg_kmers[akmer];
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
inline uint8_t DeBrujinGraph<KMERBITS>::outdegree(const bitset<SIGMA + 1>& ar) {
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
