#ifndef SCDBG_DBG_STATS_H
#define SCDBG_DBG_STATS_H

#include "dbg.h"

template<uint16_t KMERBITS>
class DeBrujinGraphStats : public DeBrujinGraph<KMERBITS> {
public:
    DeBrujinGraphStats(const uint8_t pkm, const uint32_t pc, const uint32_t psmd = 0) : DeBrujinGraph<KMERBITS>(pkm, pc,
                                                                                                                psmd) {}

    void do_stats();

private:
    string kmer_to_str(bitset<KMERBITS> kmer_str);

    void print_node(const bitset<KMERBITS>& str, uint64_t icnt, uint64_t ocnt);

    tuple<size_t, size_t, size_t>
    traverse(bitset<KMERBITS> pkmer, bool path, sparse_hash_map<bitset<KMERBITS>, uint8_t>& visited);

    tuple<size_t, size_t, size_t> traverse2(const bitset<KMERBITS>& pkmer, bool path);
};

#endif //SCDBG_DBG_STATS_H
