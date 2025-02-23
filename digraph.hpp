#pragma once
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <unordered_map>
#include <vector>

class divertex;
class diedge;

class diedge {
public:
    unsigned long int id;
    double weight = 0;
    divertex* from;
    divertex* to;
};

class divertex {
public:
    unsigned long int id;
    std::vector<diedge*> outedges;
    std::vector<diedge*> inedges;
};

class digraph {
public:
    std::vector<diedge*> edges;
    std::vector<divertex*> vertices;

    // std::unordered_map<unsigned long, unsigned long> vertex_by_id;


    unsigned long add_vertex() {
        auto* v = new divertex();
        v->id = vertices.size();
        vertices.push_back(v);
        return v->id;
    }

    unsigned long add_edge(unsigned long from, unsigned long to, double weight) {
        auto* e = new diedge();
        e->from = vertices[from];
        e->to = vertices[to];
        e->weight = weight;
        e->id = edges.size();
        edges.push_back(e);
        vertices[from]->outedges.push_back(e);
        vertices[to]->inedges.push_back(e);
        return e->id;
    }


    void print() {
        fprintf(stderr, "Graph has %zu vertices and %zu edges.\n", vertices.size(), edges.size());
    }
    ~digraph() {
        for (auto v : vertices) {
            delete v;
        }
        for (auto e : edges) {
            delete e;
        }
    }

    void bellman_ford() {
        long unsigned n = vertices.size();
        long unsigned iteration_limit = n;
        double* distances;
        long int* pred;

        distances = (double*) malloc(n * sizeof(double));
        pred = (long int*) malloc(n * sizeof(long int));

        for (long unsigned int i = 0; i < n; i++) {
            distances[i] = std::numeric_limits<double>::max();
            pred[i] = -1;
        }

        distances[0] = 0;
        pred[0] = 0;

        for (int iteration = 0; iteration < iteration_limit; iteration++) {
            fprintf(stderr, "Iteration %d.\n", iteration);
            bool update_happened = false;
            // For every edge means going through all vertices once more and listing the edges there.
           for (auto edge : edges) {
               long unsigned int from = edge->from->id;
               long unsigned int to = edge->to->id;
               double weight = edge->weight;
               if (distances[from] != std::numeric_limits<double>::max() && distances[from] + weight < distances[to]) {
                  distances[to] = distances[from] + weight;
                  pred[to] = static_cast<long int>(from);
                  update_happened = true;
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
                return;
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

        /*
        for (int i = 0; i < g.verts.size(); i++) {
            for (int j = 0; j < g.verts[i].size(); j++) {
                if (!g.verts[i][j]->reachable) {
                    continue;
                }

                for (auto edge : g.verts[i][j]->edgelist) {
                    long unsigned int from = edge->from->id;
                    long unsigned int to = edge->to->id;
                    cost_t weight = EDGE_WEIGHT(edge);
                    if (distances[from] != (cost_t)INT64_MAX && distances[from] + weight < distances[to]) {
                        // negative_cycle_found = true;
                        fprintf(stderr,
                                "Negative cycle_with_tail found in the graph. Relevant vertex with distance value %f:\n",
                                distances[from]);
                        edge->from->print(stderr);
                        fprintf(stderr, "Relevant vertex to with distance value %f:\n", distances[to]);
                        edge->to->print(stderr);
                        fprintf(stderr, "pred[from] = %ld.\n", pred[from]);

                        // Build the negative cycle_with_tail.
                        std::vector<long int> cycle_with_tail;
                        std::unordered_set<long int> visited;

                        cycle_with_tail.push_back((long int)from);
                        visited.insert((long int)from);
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
        */

        fprintf(stderr, "No negative cycles present.\n");
        free(pred);
        free(distances);
    }
};
