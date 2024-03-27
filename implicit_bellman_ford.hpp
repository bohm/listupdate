#pragma once

#include "common.hpp"
#include "permutations.hpp"
#include "implicit_graph.hpp"
#include <cstdlib>

void bellman_ford_implicit() {
    long unsigned int n = factorial(LISTSIZE) * (MEMORY::max + 1);
    fprintf(stderr, "There are %ld vertices in the graph.\n", n);

    cost_t *distances;
    long int *pred;

    distances = (cost_t *) malloc(n * sizeof(cost_t));
    pred = (long int *) malloc(n * sizeof(long int));

    for (long unsigned int i = 0; i < n; i++) {
        distances[i] = (cost_t) INT64_MAX;
        pred[i] = -1;
    }

    distances[0] = 0;
    pred[0] = 0;

    for (int iteration = 0; iteration < n; iteration++) {
        fprintf(stderr, "Iteration %d.\n", iteration);
        bool update_happened = false;
        // For every edge means going through all vertices once more and listing the edges there.
        for (long int from = 0; from < n; from++) {
            // First, recover vertex description.
            auto [v_perm, v_mem] = implicit_graph::get_vertex_information(from);
            // Go through presentation edges first.
            for (int j = 0; j < LISTSIZE; j++) {
                auto [to, weight] = implicit_graph::presentation_edge(v_perm, v_mem, j);
                if (distances[from] != (cost_t) INT64_MAX && distances[from] + weight < distances[to]) {
                    distances[to] = distances[from] + weight;
                    pred[to] = (long int) from;
                    update_happened = true;
                }
            }

            // Now repeat the BF procedure for translation edges.
            for (int j = 0; j < LISTSIZE - 1; j++) {
                auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
                if (distances[from] != (cost_t) INT64_MAX && distances[from] + weight < distances[to]) {
                    distances[to] = distances[from] + weight;
                    pred[to] = (long int) from;
                    update_happened = true;
                }
            }
        }

        if (!update_happened) {
            break;
        }
    }

    // Test for negative cycles.
    // bool negative_cycle_found = false;

    /*
    fprintf(stderr, "[");
    for (long int x = 0; x < n; x++) {
        fprintf(stderr, "%ld,", distances[x]);
    }
         fprintf(stderr, "]\n");

     */

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
    free(pred);
    free(distances);
}
