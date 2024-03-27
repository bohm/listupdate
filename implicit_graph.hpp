#pragma once
#include <vector>

#include "common.hpp"
#include "algorithm.hpp"
#include "permutations.hpp"

class implicit_graph {
public:

    static std::pair<permutation, memory_pairs> get_vertex_information(long int vertex) {
        uint64_t memory_section = vertex % (MEMORY::max + 1);
        uint64_t permutation_section = vertex / (MEMORY::max + 1);
        memory_pairs ret2; ret2.data = memory_section;
        permutation ret1 = perm_from_index_quadratic(permutation_section);
        return {ret1, ret2};
    }

    static std::pair<long int, cost_t> presentation_edge(permutation v_perm, memory_pairs v_mem, int presented_item ) {

        int alg_cost = ALG_SINGLE_STEP(&v_perm, &v_mem, presented_item); // v_perm will get edited.
        int opt_cost = presented_item;
        long int target = lexindex_quadratic(&v_perm) * (MEMORY::max + 1) + v_mem.data;
        return {target, EDGE_WEIGHT(opt_cost, alg_cost)};
    }

    static std::pair<long int, cost_t> translation_edge(permutation v_perm, memory_pairs v_mem,
                                                  int translation_index) {
        permutation single_swap = IDENTITY;
        swap(&single_swap, translation_index);

        MEMORY mem_copy = v_mem.recompute(&single_swap);
        recompute_alg_perm(&v_perm, &single_swap);
        long int target = lexindex_quadratic(&v_perm) * (MEMORY::max + 1) + mem_copy.data;
        return {target, EDGE_WEIGHT(1, 0)};
    }
};