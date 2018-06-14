#ifndef SCDBG_COMPRESSED_DBG_H
#define SCDBG_COMPRESSED_DBG_H

#include <iostream>
#include <fstream>
#include <ostream>
#include <array>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>

#include "config.h"
#include "utils.hpp"

class CompressedDeBrujinGraph {
public:
    CompressedDeBrujinGraph(const std::string& fname);

private:
    bool load_dbg(std::ifstream& f);

    bool load_edge_list(std::ifstream& f);

    bool load_BL_BF(std::ifstream& f, size_t length, sdsl::rank_support_v<>& rs,
                    sdsl::select_support_mcl<>& ss);

    void load_bit_vectors(std::ifstream& f);

    size_t forward(size_t index);

    size_t num_of_edges = 0;
    size_t num_of_nodes = 0;
    std::array<size_t, SIGMA> F = {};
    size_t num_of_color_classes = 0;
    // sdsl::int_vector<8> edge_list;
    sdsl::wt_blcd<> edge_list;
    // sdsl::bit_vector BL;
    sdsl::rank_support_v<> BL_rs;
    sdsl::select_support_rrr<> BL_ss;
    // sdsl::bit_vector BF;
    sdsl::rank_support_v<> BF_rs;
    sdsl::select_support_mcl<> BF_ss;
    std::ifstream color_classes_file;
};


#endif //SCDBG_COMPRESSED_DBG_H
