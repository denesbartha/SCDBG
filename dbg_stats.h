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

    void print_node(const bitset<KMERBITS> &str, uint64_t icnt, uint64_t ocnt);

    Roaring get_color(const bitset<KMERBITS> &pkmer);
};

#endif //SCDBG_DBG_STATS_H
