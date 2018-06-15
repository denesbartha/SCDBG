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
    class BV {
    public:
        // BV(const sdsl::bit_vector& bv) : rs(&bv), ss(&bv) {}
        //
        // BV(const BV& bv) : rs(bv.rs), ss(bv.ss) {}

        void resize(const size_t length) { bv.resize(length); }

        sdsl::int_vector<1>::reference operator[](const size_t i) { return bv[i]; }

        void init_support(const std::array<size_t, SIGMA>& F) {
            sdsl::util::init_support(rs, &bv);
            sdsl::util::init_support(ss, &bv);

            for (auto i = 0; i < F.size(); ++i) {
                F_rank[i] = rs.rank(F[i] + 1);
            }
        }

        size_t rank(size_t i) const { return rs.rank(i); }

        size_t select(size_t i) const { return ss.select(i); }

        size_t size() const { return rs.size(); }

    private:
        sdsl::bit_vector bv;
        sdsl::rank_support_v<> rs;
        sdsl::select_support_mcl<> ss;
        // culmulative sums of rank values for table F
        std::array<size_t, SIGMA> F_rank = { };
    } BL, BF;

    bool load_dbg(std::ifstream& f);

    bool load_edge_list(std::ifstream& f);

    bool load_BL_BF(std::ifstream& f, size_t length, BV& bv);

    void load_bit_vectors(std::ifstream& f);

    size_t forward(size_t index);

    size_t num_of_edges = 0;
    size_t num_of_nodes = 0;
    std::array<size_t, SIGMA> F = {};
    size_t num_of_color_classes = 0;
    sdsl::wt_blcd<> edge_list;
    std::ifstream color_classes_file;
};

#endif //SCDBG_COMPRESSED_DBG_H
