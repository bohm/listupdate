#pragma once

#include <cstdlib>
#include <omp.h>
#include <atomic>
#include "common.hpp"
#include "permutations.hpp"
#include "implicit_graph.hpp"
#include "storage.hpp"

void bellman_ford_compare_exchange() {
    long unsigned int n = factorial(LISTSIZE) * (MEMORY::max + 1);
    fprintf(stderr, "There are %ld vertices in the graph.\n", n);

    std::atomic<cost_t> *distances;
    distances = new std::atomic<cost_t>[n];
    // The predecessor code is not implemented for the implicit version. This is unfortunate,
    // but gives out large memory savings.

#pragma omp parallel for
    for (long unsigned int i = 0; i < n; i++) {
        distances[i] = (cost_t) INT32_MAX;
    }

    distances[0] = 0;

    for (int iteration = 0; iteration < n; iteration++) {
        fprintf(stderr, "Iteration %d.\n", iteration);
        bool update_happened = false;

        fprintf(stderr, "Iteration %d distances[0]: %f.\n", iteration, distances[0].load());


#pragma omp parallel for
        // For every edge means going through all vertices once more and listing the edges there.
        for (long int from = 0; from < n; from++) {
            // First, recover vertex description.
            auto [v_perm, v_mem] = implicit_graph::get_vertex_information(from);
            // Go through presentation edges first.
            for (int j = 0; j < LISTSIZE; j++) {
                auto [to, weight] = implicit_graph::presentation_edge(v_perm, v_mem, j);
                // We do not really care about the from number, which could update on the fly. The problem
                // is that distances[to] should not change while we do the comparison.

                cost_t dist_from = distances[from];
                  if (dist_from != (cost_t) INT32_MAX) {
                    cost_t dist_to = distances[to].load(std::memory_order_relaxed);
                    if (dist_from + weight < dist_to) {
                        update_happened = true;
                    }

                    while (dist_from + weight < dist_to) {
                        distances[to].compare_exchange_weak(dist_to, dist_from + weight);
                        // dist_to = distances[to].load(std::memory_order_relaxed);
                    }
                }
            }

            // Now repeat the BF procedure for translation edges.
            for (int j = 0; j < LISTSIZE - 1; j++) {
                auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
                cost_t dist_from = distances[from];
                if (dist_from != (cost_t) INT32_MAX) {
                    cost_t dist_to = distances[to].load(std::memory_order_relaxed);
                    if (dist_from + weight < dist_to) {
                        update_happened = true;
                    }

                    while (dist_from + weight < dist_to) {
                        distances[to].compare_exchange_weak(dist_to, dist_from + weight);
                        // dist_to = distances[to].load(std::memory_order_relaxed);
                    }
                }
            }
        }

        if (!update_happened) {
            break;
        }

        if (distances[0] < 0.0) {
            fprintf(stderr, "Negative cycle found in the graph.\n");
            // print_array(n, reinterpret_cast<cost_t *>(distances));
            // write_distance_array(reinterpret_cast<cost_t *>(distances), n);
            // implicit_graph::print_negative_cycle(n, reinterpret_cast<cost_t *>(distances));
            return;
        }
    }

    // Test for negative cycles.
    // bool negative_cycle_found = false;

    // fprintf(stderr, "Iteration %d.\n", iteration);
    bool update_happened = false;
    // For every edge means going through all vertices once more and listing the edges there.
    for (long int from = 0; from < n; from++) {
        // First, recover vertex description.
        auto [v_perm, v_mem] = implicit_graph::get_vertex_information(from);
        // Go through presentation edges first.
        for (int j = 0; j < LISTSIZE; j++) {
            auto [to, weight] = implicit_graph::presentation_edge(v_perm, v_mem, j);
            if (distances[from] != (cost_t) INT64_MAX && distances[from] + weight < distances[to]) {
                fprintf(stderr, "Negative cycle found in the graph.\n");
                return;
            }
        }

        // Now repeat the BF procedure for translation edges.
        for (int j = 0; j < LISTSIZE - 1; j++) {
            auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
            if (distances[from] != (cost_t) INT64_MAX && distances[from] + weight < distances[to]) {
                fprintf(stderr, "Negative cycle found in the graph.\n");
                return;
            }
        }
    }

    fprintf(stderr, "No negative cycles present.\n");
    // print_array(n, reinterpret_cast<cost_t *>(distances));
    // write_distance_array(reinterpret_cast<cost_t *>(distances), n);

    delete distances;
}
