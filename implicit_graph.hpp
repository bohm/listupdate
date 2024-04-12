#pragma once
#include <vector>
#include <unordered_set>
#include <algorithm>

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

    static bool presentation_predecessor(long int parent_candidate, long int child_candidate) {
        auto [v_perm, v_mem] = get_vertex_information(parent_candidate);
        for (int j = 0; j < LISTSIZE; j++) {
            auto [to, _] = implicit_graph::presentation_edge(v_perm, v_mem, j);
            if (to == child_candidate) {
                return true;
            }
        }
        return false;
    }


    static bool translation_predecessor(long int parent_candidate, long int child_candidate) {
        auto [v_perm, v_mem] = get_vertex_information(parent_candidate);
        for (int j = 0; j < LISTSIZE - 1; j++) {
            auto [to, _] = implicit_graph::translation_edge(v_perm, v_mem, j);
            if (to == child_candidate) {
                return true;
            }
        }
        return false;
    }

    static long int linear_time_predecessor(long int distances_size, cost_t *distances, long int vertex) {
        cost_t distance_v = distances[vertex];
        for (long int parent_cand = 0; parent_cand < distances_size; parent_cand++) {
            if (parent_cand == vertex) {
                continue;
            }
            auto [v_perm, v_mem] = get_vertex_information(parent_cand);
            // Presentation edge check.
            for (int j = 0; j < LISTSIZE; j++) {
                auto [to, weight] = implicit_graph::presentation_edge(v_perm, v_mem, j);
                if (to == vertex) {
                    if (distance_v == distances[parent_cand] + weight) {
                        return parent_cand;
                    }
                }
            }

            // Translation edge check.
            for (int j = 0; j < LISTSIZE - 1; j++) {
                auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
                if (to == vertex) {
                    if (distance_v == distances[parent_cand] + weight) {
                        return parent_cand;
                    }
                }
            }
        }

        return -1;
    }

    static void vertex_print(long int id, permutation *p, memory_pairs *mem, FILE *f) {

        fprintf(f, "%ld [label=\"%lu,", id, mem->data);
        print_permutation(p, f, false);
        fprintf(f, "\"];\n");

    }

    static void locate_edge_and_print(long int parent, long int child) {

        auto [v_perm, v_mem] = get_vertex_information(parent);
        // Presentation edge check.
        for (int j = 0; j < LISTSIZE; j++) {
            auto [to, weight] = implicit_graph::presentation_edge(v_perm, v_mem, j);
            if (to == child) {
                fprintf(stderr, "%lu -> %lu [label=\"req: %d, edge_weight %f\"];\n",
                        parent, child, j, weight);
                return;
            }
        }

        // Translation edge check.
        for (int j = 0; j < LISTSIZE - 1; j++) {
            auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
            if (to == child) {
                fprintf(stderr, "%lu -> %lu [label=\"swap %d,%d\"];\n", parent, child, j, j+1);
                return;
            }
        }

        fprintf(stderr, "Error: edge between %ld and %ld not located.\n", parent, child);
        exit(-1);
    }


    static void print_vertex_sequence(std::vector<long int> seq) {
        for (int counter = 0; counter < seq.size(); counter++) {
            fprintf(stderr, "Vertex %d/%zu:\n", counter, seq.size());
            auto [v_perm, v_mem] = get_vertex_information(seq[counter]);

            vertex_print(seq[counter], &v_perm, &v_mem, stderr);
            fprintf(stderr, "Memory content for vertex %d/%zu: ", counter, seq.size());
            v_mem.full_print();
            fprintf(stderr, "\n");
            if (counter < seq.size() - 1) {
                locate_edge_and_print(seq[counter], seq[counter+1]);
            }
        }
    }

    static void print_negative_cycle(long int distances_size, cost_t* distances, long int from = 0) {
        // Build the negative cycle.
        std::vector<long int> cycle;
        std::unordered_set<long int> visited;

        cycle.push_back((long int) from);
        visited.insert((long int) from);
        long int p = linear_time_predecessor(distances_size, distances, from);
        while (!visited.contains(p)) {
            cycle.push_back(p);
            visited.insert(p);
            p = linear_time_predecessor(distances_size, distances, p);
        }
        cycle.push_back(p);
        fprintf(stderr, "One negative sequence (cycle with tail) has length %zu.\n", cycle.size());
        std::reverse(cycle.begin(), cycle.end());
        print_vertex_sequence(cycle);
    }
};