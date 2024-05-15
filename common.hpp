#pragma once
#include <array>
#include <cstdint>
#include <cstdio>

// using cost_t = long int;
// using cost_t = long double; // Makes the most sense, to get the best precision if memory is not an issue.
using cost_t = float; // If memory is an issue, you can use this.

constexpr unsigned short LISTSIZE = 4;

#define MEMORY memory_pairs
#define ALG_SINGLE_STEP alg_single_step_xoror

// #define MEMORY memory_bitfield
// #define ALG_SINGLE_STEP alg_single_step_bitfield

constexpr float EPSILON = 0.0001;
// constexpr long double RATIO = 3.6667;
constexpr long double RATIO = 3.15;
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


void print_array(unsigned long int len, cost_t *array) {
    fprintf(stderr, "[");
    for (long int x = 0; x < len; x++) {
        fprintf(stderr, "%f,", array[x]);
    }
    fprintf(stderr, "]\n");
}