#pragma once
#include <array>
#include <cstdint>
#include <cstdio>

#include "constants.hpp"
// using cost_t = long int;
// using cost_t = long double; // Makes the most sense, to get the best precision if memory is not an issue.
using cost_t = float; // If memory is an issue, you can use this.



// #define MEMORY memory_pairs
// #define ALG_SINGLE_STEP alg_single_step_xoror

//#define MEMORY memory_perm

//#define ALG_SINGLE_STEP alg_single_step_mru_first_inversion
// #define ALG_SINGLE_STEP alg_single_step_lessrecent
// #define ALG_INFO alg_single_step_mru_eager_info
#define MEMORY memory_bitfield
#define ALG_SINGLE_STEP alg_single_step_bitfield

constexpr bool ALG_DEBUG = false;
constexpr bool GRAPH_DEBUG = false;
constexpr bool FRONT_ACCESS_COSTS_ONE = true;


#define EDGE_WEIGHT edge_weight_param


cost_t edge_weight_param(cost_t opt_cost, cost_t alg_cost) {
    return ((cost_t) RATIO)*opt_cost - alg_cost;
}

using array_as_permutation = std::array<short, LISTSIZE>;

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


void print_array(unsigned long int len, short *array) {
    fprintf(stderr, "[");
    for (long int x = 0; x < len; x++) {
        fprintf(stderr, "%hu,", array[x]);
    }
    fprintf(stderr, "]\n");
}

constexpr uint64_t fact(uint64_t n) {
    return n <= 1 ? 1 : (n* fact(n-1));
}

static constexpr std::array<uint64_t, TESTSIZE+1> fact_array() {
    std::array<uint64_t, TESTSIZE+1> ret{};
    for (int i = 0; i <= TESTSIZE; i++) {
        ret[i] = fact(i);
    }
    return ret;
}

constexpr std::array<uint64_t, TESTSIZE+1> factorial = fact_array();

constexpr short diameter_bound(short n) {
    return (n*(n-1))/2;
}

inline bool triple_contains(const std::array<int8_t, 3>* ar, const short el) {
    return (el == (*ar)[0] || el == (*ar)[1] || el == (*ar)[2]);
}