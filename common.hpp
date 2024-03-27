#pragma once
#include <array>

// using cost_t = long int;
using cost_t = long double;

constexpr unsigned short LISTSIZE = 5;

#define MEMORY memory_pairs
#define ALG_SINGLE_STEP alg_single_step_original

// #define MEMORY memory_bitfield
// #define ALG_SINGLE_STEP alg_single_step_bitfield


constexpr long double RATIO = 3.6667;
#define EDGE_WEIGHT edge_weight_param


cost_t edge_weight_param(cost_t opt_cost, cost_t alg_cost) {
    return ((cost_t) RATIO)*opt_cost - alg_cost;
}

using permutation = std::array<short, LISTSIZE>;

constexpr std::array<int, LISTSIZE*LISTSIZE> canonical_ordering() {
    std::array<int, LISTSIZE*LISTSIZE> co {0};
    int counter = 0;
    for (int i = 0; i <LISTSIZE; i++) {
        for (int j = i+1; j < LISTSIZE; j++) {
            co[i*LISTSIZE+j] = counter++;
        }
    }
    return co;
}

constexpr std::array<int, LISTSIZE*LISTSIZE> canonical_order = canonical_ordering();

constexpr uint64_t max_memory_pairs = (1LLU << (canonical_order[(LISTSIZE - 2) * LISTSIZE + (LISTSIZE - 1)] + 1)) - 1;