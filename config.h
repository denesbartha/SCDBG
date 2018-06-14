#ifndef SCDBG_CONFIG_H
#define SCDBG_CONFIG_H

// SIGMA (for DNA, it is 4...)
#define SIGMA       4u

//log_2(SIGMA + 1)
#define LOGSIGMA    3u

// maximum number of colors
#define MAXCOLORS   100u // 100u 569u 95146u

static const char base[5] = {'$', 'A', 'C', 'G', 'T'};

#endif //SCDBG_CONFIG_H
