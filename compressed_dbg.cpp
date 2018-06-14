#include "compressed_dbg.h"

using namespace std;
using namespace sdsl;

// constructor for loading DBG file and colour information
CompressedDeBrujinGraph::CompressedDeBrujinGraph(const std::string& fname) {
    // .dbg
    ifstream dbg_f(fname + ".dbg", ios::in | ios::binary);
    if (dbg_f.is_open()) {
        cerr << "Loading DBG file " << fname << ".dbg..." << endl;
        if (!load_dbg(dbg_f)) {
            return;
        }
    }
    else {
        cerr << "Unable to open file " << fname << ".dbg..." << endl;
        return;
    }
    dbg_f.close();

    // .color_classes - don't load the whole file into the memory (it could be really large...), just store file handler
    color_classes_file.open(fname + ".color_classes", ios::in | ios::binary);
    if (color_classes_file.is_open()) {
        // read the number of color classes
        cerr << "Loading color classes file " << fname << ".color_classes..." << endl;
        color_classes_file.read(reinterpret_cast<char *>(&num_of_color_classes), sizeof(num_of_color_classes));
    }
    else {
        cerr << "Unable to open file " << fname << ".color_classes..." << endl;
        return;
    }

    // .bit_vectors
    ifstream bit_vectors_f(fname + ".bit_vectors", ios::in | ios::binary);
    if (bit_vectors_f.is_open()) {
        cerr << "Loading bit vectors file " << fname << ".bit_vectors..." << endl;
        load_bit_vectors(bit_vectors_f);
    }
    else {
        cerr << "Unable to open file " << fname << ".bit_vectors..." << endl;
        return;
    }
    bit_vectors_f.close();
}


// loads the given DBG file - labels, F, B_L, B_F
bool CompressedDeBrujinGraph::load_dbg(std::ifstream& f) {
    // read the number of edges
    f.read(reinterpret_cast<char *>(&num_of_edges), sizeof(num_of_edges));
    // load the labels
    if (!load_edge_list(f)) {
        return false;
    }

    // load table F
    cerr << "Loading table F..." << endl;
    for (auto i = 0; i < SIGMA; ++i) {
        f.read(reinterpret_cast<char *>(&F[i]), sizeof(F[i]));
    }

    // load B_L
    cerr << "Loading bit vector B_L..." << endl;
    load_BL_BF(f, num_of_edges, BL_rs, BF_ss);
    // load B_F - note that it's size is (# of edges) - 1
    cerr << "Loading bit vector B_F..." << endl;
    load_BL_BF(f, num_of_edges - 1, BF_rs, BF_ss);

    return true;
}


bool CompressedDeBrujinGraph::load_edge_list(std::ifstream& f) {
    char *buffer = nullptr;
    cerr << "Loading edge list..." << endl;
    load_from_file(f, buffer, LOGSIGMA * num_of_edges);
    if (buffer == nullptr) {
        cerr << "Error loading edge list...";
        return false;
    }

    cerr << "Processing edge list..." << endl;
    // put the edges into a wavelet tree
    // edge_list.resize(num_of_edges);
    size_t buffer_index = 0;
    // for (size_t i = 0; i < num_of_edges; ++i) {
    //     unsigned char ac = 0;
    //     for (uint8_t j = 0; j < LOGSIGMA; ++j, ++buffer_index) {
    //         ac |= buffer[buffer_index / 8] & (1 << (j + (buffer_index % 8)));
    //     }
    //     edge_list[i] = ac;
    // }
    // construct(edge_list, )
    // reset_buffer(buffer);
    // rank_support_il<1, 8>(edge_list);

    // sdsl::rrr_vector<3> rr;
    // rank_support_rrr rs(rr);

    return true;
}


bool CompressedDeBrujinGraph::load_BL_BF(std::ifstream& f, size_t length, sdsl::rank_support_rrr<>& rank_support,
                                         sdsl::select_support_rrr<>& select_support_rrr) {
    char *buffer = nullptr;
    load_from_file(f, buffer, length);
    if (buffer == nullptr) {
        cerr << "Error loading bit vector...";
        return false;
    }

    cerr << "Processing edge list..." << endl;
    // put the edges into a wavelet tree
    bit_vector bv(length, 0);
    size_t buffer_index = 0;
    for (size_t i = 0; i < num_of_edges; ++i, ++buffer_index) {
        bv[i] = (buffer[buffer_index / 8] & (1 << (buffer_index % 8))) != 0;
    }
    reset_buffer(buffer);
    // rank_support_v<> rs(&bv);
    // rank_support = rs;

    return true;
}


// loads the bit vectors - store, label, boundary
void CompressedDeBrujinGraph::load_bit_vectors(ifstream& f) {

}


size_t CompressedDeBrujinGraph::forward(size_t index) {
    // get the current character from L
    // get the rank of the character c
    // find the next 1 in B_F - get the rank of this (extract the rank of the first c)
    // select the above calculated value in B_L => find the node's first position by reading 0s backward
}

