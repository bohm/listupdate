#pragma once

#include "common.hpp"
#include "permutations.hpp"
#include "implicit_graph.hpp"
#include <cstdlib>
#include <omp.h>
#include <atomic>

void bellman_ford_implicit() {
    long unsigned int n = factorial(LISTSIZE) * (MEMORY::max + 1);
    fprintf(stderr, "There are %ld vertices in the graph.\n", n);

    cost_t *distances;

    distances = (cost_t *) malloc(n * sizeof(cost_t));
    // The predecessor code is not implemented for the implicit version. This is unfortunate,
    // but gives out large memory savings.

    for (long unsigned int i = 0; i < n; i++) {
        distances[i] = (cost_t) INT32_MAX;
    }

    distances[0] = 0;

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
                if (distances[from] != (cost_t) INT32_MAX && distances[from] + weight < distances[to]) {
                    distances[to] = distances[from] + weight;
                    update_happened = true;
                }
            }

            // Now repeat the BF procedure for translation edges.
            for (int j = 0; j < LISTSIZE - 1; j++) {
                auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
                if (distances[from] != (cost_t) INT32_MAX && distances[from] + weight < distances[to]) {
                    distances[to] = distances[from] + weight;
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
    free(distances);
}


cost_t sum_of_negative_edges() {
        long unsigned int n = factorial(LISTSIZE) * (MEMORY::max + 1);
	std::atomic<float> sum = 0.0;
        // For every edge means going through all vertices once more and listing the edges there.
        #pragma omp parallel for
        for (long int from = 0; from < n; from++) {
            // First, recover vertex description.
            auto [v_perm, v_mem] = implicit_graph::get_vertex_information(from);
            // Go through presentation edges first.
            for (int j = 0; j < LISTSIZE; j++) {
                auto [to, weight] = implicit_graph::presentation_edge(v_perm, v_mem, j);
		if (weight < 0.0) {
			sum += weight;
		}
	    }
    
            for (int j = 0; j < LISTSIZE - 1; j++) {
                auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
		if (weight < 0.0) {
			sum += weight;
		}
            }
	}

	fprintf(stderr, "The total sum on all negative edges is %f.\n", sum.load());
	return  (cost_t) sum;	
}

void bellman_ford_implicit_parallel() {
    long unsigned int n = factorial(LISTSIZE) * (MEMORY::max + 1);
    fprintf(stderr, "There are %ld vertices in the graph.\n", n);

    cost_t *distances;
    omp_lock_t locks[1024*1024];
    for (int i = 0; i <1024*1024; i++) {
        omp_init_lock(&(locks[i]));
    }

    distances = (cost_t *) malloc(n * sizeof(cost_t));
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

	// This does not belong to Bellman-Ford, and it should be commented out when running
	// it in performance mode.
	cost_t min_distance = INT32_MAX;

	for (long int i = 0; i < n; i++) {
        	min_distance = std::min(min_distance, distances[i]);
    	}

	fprintf(stderr, "Iteration %d min distance: %f.\n", iteration, min_distance);
	fprintf(stderr, "Iteration %d distances[0]: %f.\n", iteration, distances[0]);

 
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
                long int lock_num = to % (1024*1024);
                omp_set_lock(&locks[lock_num]);
                cost_t dist_to = distances[to];
                if (dist_from != (cost_t) INT32_MAX && dist_from + weight < dist_to) {
                    distances[to] = dist_from + weight;
                    update_happened = true;
                }
                omp_unset_lock(&locks[lock_num]);
            }

            // Now repeat the BF procedure for translation edges.
            for (int j = 0; j < LISTSIZE - 1; j++) {
                auto [to, weight] = implicit_graph::translation_edge(v_perm, v_mem, j);
                cost_t dist_from = distances[from];
                long int lock_num = to % (1024*1024);
                omp_set_lock(&locks[lock_num]);
                cost_t dist_to = distances[to];
                if (dist_from != (cost_t) INT32_MAX && dist_from + weight < dist_to) {
                    distances[to] = dist_from + weight;
                    update_happened = true;
                }
                omp_unset_lock(&locks[lock_num]);
            }
        }

        if (!update_happened) {
            break;
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
    free(distances);
}
