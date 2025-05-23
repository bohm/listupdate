#pragma once

#include <cstdlib>
#include "common.hpp"
#include "old_perm_functions.hpp"
#include "graph.hpp"


// Note: Bellman for now has a reachability clause. This only makes sense if reachability is computed
// and would slow us down if everything is reachable.

void bellman_ford() {
    long unsigned int n = factorial[LISTSIZE] * (MEMORY::max + 1);
    fprintf(stderr, "There are %ld vertices in the graph.\n", n);
    // long unsigned iteration_limit = n;
    long unsigned iteration_limit = g.reachable_vertices;
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

    for (int iteration = 0; iteration < iteration_limit; iteration++) {
        fprintf(stderr, "Iteration %d.\n", iteration);
        bool update_happened = false;
        // For every edge means going through all vertices once more and listing the edges there.
        for (int i = 0; i < g.verts.size(); i++) {
            for (int j = 0; j < g.verts[i].size(); j++) {
                if (!g.verts[i][j]->reachable) {
                    continue;
                }

                for (auto edge: g.verts[i][j]->edgelist) {
                    long unsigned int from = edge->from->id;
                    long unsigned int to = edge->to->id;
                    cost_t weight = EDGE_WEIGHT(edge);
                    if (distances[from] != (cost_t) INT64_MAX && distances[from] + weight < distances[to]) {
                        distances[to] = distances[from] + weight;
                        pred[to] = (long int) from;
                        update_happened = true;
                    }
                }
            }
        }

        if (!update_happened) {
            fprintf(stderr, "No negative cycles present.\n");
            free(pred);
            free(distances);
            return;
        }

        if (distances[0] < 0.0) {
            fprintf(stderr, "Negative cycle found in the graph.\n");
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

    for (int i = 0; i < g.verts.size(); i++) {
        for (int j = 0; j < g.verts[i].size(); j++) {
            if (!g.verts[i][j]->reachable) {
                continue;
            }

            for (auto edge: g.verts[i][j]->edgelist) {
                long unsigned int from = edge->from->id;
                long unsigned int to = edge->to->id;
                cost_t weight = EDGE_WEIGHT(edge);
                if (distances[from] != (cost_t) INT64_MAX && distances[from] + weight < distances[to]) {
                    // negative_cycle_found = true;
                    fprintf(stderr, "Negative cycle_with_tail found in the graph. Relevant vertex with distance value %f:\n",
                            distances[from]);
                    edge->from->print(stderr);
                    fprintf(stderr, "Relevant vertex to with distance value %f:\n", distances[to]);
                    edge->to->print(stderr);
                    fprintf(stderr, "pred[from] = %ld.\n", pred[from]);

                    // Build the negative cycle_with_tail.
                    std::vector<long int> cycle_with_tail;
                    std::unordered_set<long int> visited;

                    cycle_with_tail.push_back((long int) from);
                    visited.insert((long int) from);
                    long int p = pred[from];
                    while (!visited.contains(p)) {
                        cycle_with_tail.push_back(p);
                        visited.insert(p);
                        p = pred[p];
                    }
                    cycle_with_tail.push_back(p);
                    std::reverse(cycle_with_tail.begin(), cycle_with_tail.end());
                    // Until here, we get a cycle with a tail. We clean off the tail to have just the cycle.
                    std::vector<long int> cycle;
                    cycle.push_back(cycle_with_tail[0]);
                    long int i = 1;
                    while (cycle_with_tail[i] != cycle[0]) {
                        cycle.push_back(cycle_with_tail[i++]);
                    }
                    cycle.push_back(cycle[0]);
                    fprintf(stderr, "One negative sequence (cycle) has length %zu.\n", cycle.size());
                    double acost = total_alg_cost(cycle);
                    double ocost = total_opt_cost(cycle);
                    fprintf(stderr, "alg cost: %F, opt cost %F, ratio %F.\n", acost, ocost, acost / ocost);
                    print_vertex_sequence(cycle);

                    free(pred);
                    free(distances);
                    return;
                }
            }
        }
    }
    fprintf(stderr, "No negative cycles present.\n");
    free(pred);
    free(distances);
}
