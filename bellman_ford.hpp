#pragma once

#include <cstdlib>
#include "common.hpp"
#include "permutations.hpp"
#include "graph.hpp"

void bellman_ford() {
    long unsigned int n = factorial(LISTSIZE)*(MEMORY::max + 1);
    fprintf(stderr, "There are %ld vertices in the graph.\n", n);

    cost_t *distances;
    long int *pred;

    distances = (cost_t*) malloc(n*sizeof(cost_t));
    pred = (long int*) malloc(n*sizeof(long int));

    for (long unsigned int i = 0; i < n; i ++) {
        distances[i] = (cost_t) INT64_MAX;
        pred[i] = -1;
    }

    distances[0] = 0;
    pred[0] = 0;

    for (int iteration = 0; iteration < n; iteration++) {
        // fprintf(stderr, "Iteration %d.\n", iteration);
        bool update_happened = false;
        // For every edge means going through all vertices once more and listing the edges there.
        for (int i = 0; i < g.verts.size(); i++) {
            for (int j = 0; j < g.verts[i].size(); j++) {
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
            for (auto edge: g.verts[i][j]->edgelist) {
                long unsigned int from = edge->from->id;
                long unsigned int to = edge->to->id;
                cost_t weight = EDGE_WEIGHT(edge);
                if (distances[from] != (cost_t) INT64_MAX && distances[from] + weight < distances[to]) {
                    // negative_cycle_found = true;
                    fprintf(stderr, "Negative cycle found in the graph. Relevant vertex with distance value %Lf:\n",
                            distances[from]);
                    edge->from->print(stderr);
                    fprintf(stderr, "Relevant vertex to with distance value %Lf:\n", distances[to]);
                    edge->to->print(stderr);
                    fprintf(stderr, "pred[from] = %ld.\n", pred[from]);

                    // Build the negative cycle.
                    std::vector<long int> cycle;
                    std::unordered_set<long int> visited;

                    cycle.push_back((long int) from);
                    visited.insert((long int) from);
                    long int p = pred[from];
                    while(!visited.contains(p)) {
                        cycle.push_back(p);
                        visited.insert(p);
                        p = pred[p];
                    }
                    cycle.push_back(p);
                    fprintf(stderr, "One negative sequence (cycle with tail) has length %zu.\n", cycle.size());
                    std::reverse(cycle.begin(), cycle.end());
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
