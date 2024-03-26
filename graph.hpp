#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <limits>
#include <unordered_set>
#include "memory.hpp"
#include "permutations.hpp"
#include "algorithm.hpp"

class adversary_vertex;

using vert_container = std::array<adversary_vertex*, max_memory+1>;

class graph;

class graph {
public:
    uint64_t edgecounter = 0;
    std::array<vert_container, factorial(LISTSIZE)> verts;

    adversary_vertex *get_vert(permutation *perm, memory m) const
    {
        return verts[lexindex_quadratic(perm)][m.data];
    }

    adversary_vertex *get_vert(long int id);

};

graph g;

class adv_outedge {
public:
    adversary_vertex *from = nullptr;
    adversary_vertex *to = nullptr;
    short presented_item = -1;
    cost_t alg_cost = 0;
    cost_t opt_cost = 0;
    int opt_swap = 0; // encode a single swap by the target position, so 0...LISTSIZE-2.
    uint64_t id = 0;

    adv_outedge(uint64_t given_id, adversary_vertex* f, adversary_vertex *t, short p_item,
                cost_t a_cost, cost_t o_cost) {
        id = given_id;
        from = f;
        to = t;
        presented_item = p_item;
        alg_cost = a_cost;
        opt_cost = o_cost;
        opt_swap = -1;
    }

    // Translation edge.
    adv_outedge(uint64_t given_id, adversary_vertex *f, adversary_vertex *t, int o_swap) {
        id = given_id;
        from = f;
        to = t;
        alg_cost = 0;
        opt_cost = 1;
        opt_swap = o_swap;
    }

    void print(FILE *f);
};

// A vertex before OPT presents an item.
class adversary_vertex {
public:
    uint64_t id;
    permutation perm;
    memory mem;
    std::vector<adv_outedge*> edgelist = {};


    adversary_vertex(permutation *p, memory m) {
        perm = *p;
        mem = m;
        id = lexindex_quadratic(&perm) * (max_memory+1) + mem.data;
    }

    inline int position(short item) const {
        for (int i = 0; i < LISTSIZE; i++) {
            if (perm[i] == item) {
                return i;
            }
        }
        return 0;
    }

    void build_presentation_edges() {
        for (short item = 0; item < LISTSIZE; item++) {
            permutation perm_copy(perm);
            memory mem_copy(mem);
            int alg_cost = ALG_SINGLE_STEP(&perm_copy, &mem_copy, item);
            int opt_cost = item;
            adversary_vertex *target = g.get_vert(&perm_copy, mem_copy);
            auto *edge = new adv_outedge(g.edgecounter++, this, target, item, alg_cost, opt_cost);
            edgelist.push_back(edge);
        }
    }

    void build_translation_edges() {
        for (int opt_swap = 0; opt_swap < LISTSIZE-1; opt_swap++) {
            permutation single_swap = IDENTITY;
            swap(&single_swap, opt_swap);

            permutation perm_copy(perm);
            memory mem_copy = recompute_memory(mem, &single_swap);
            recompute_alg_perm(&perm_copy, &single_swap);
            adversary_vertex *target = g.get_vert(&perm_copy, mem_copy);
            auto *edge = new adv_outedge(g.edgecounter++, this, target, opt_swap);
            edgelist.push_back(edge);
        }
    }


    void print(FILE *f) {

        fprintf(f, "%lu [label=\"%lu,", id, mem.data);
        print_permutation(&perm, f, false);
        fprintf(f, "\";\n");

    }

    void print_edges(FILE *f) {
        for (auto e: edgelist) {
            e->print(f);
        }
    }


};


void adv_outedge::print(FILE* f) {
    if (presented_item == -1) {
        fprintf(f, "%lu -> %lu [label=\"swap %d,%d\";\n", from->id, to->id, opt_swap, opt_swap+1);
    } else {
        fprintf(f, "%lu -> %lu [label=\"req: %d, a_cost: %Lf, o_cost: %Lf\"];\n", from->id, to->id, presented_item,
                alg_cost, opt_cost);
    }

}

void add_vertex_to_graph(permutation *perm, memory m) {
    auto *v = new adversary_vertex(perm, m);
    g.verts[lexindex_quadratic(perm)][m.data] = v;
}


adversary_vertex* graph::get_vert(long int id) {
        uint64_t memory_section = id % (max_memory+1);
        uint64_t permutation_section = id / (max_memory+1);
        assert((long int) verts[permutation_section][memory_section]->id == id);
        return verts[permutation_section][memory_section];
}

void create_graph() {
    // Build vertices.
    iterate_over_memory_and_permutation(add_vertex_to_graph);

    // Build edges. There are two types: presentation edges and translation edges (OPT swaps).
    for (int i = 0; i < g.verts.size(); i++) {
        for (int j = 0; j < g.verts[i].size(); j++) {
            g.verts[i][j]->build_presentation_edges();
            g.verts[i][j]->build_translation_edges();
        }
    }
};

void print_graph(FILE *f) {
    for (int i = 0; i < g.verts.size(); i++) {
        for (int j = 0; j < g.verts[i].size(); j++) {
            g.verts[i][j]->print(f);
            g.verts[i][j]->print_edges(f);
        }
    }
}

constexpr long double RATIO = 3.6667;
#define EDGE_WEIGHT edge_weight_param

cost_t edge_weight_param(adv_outedge *e) {
    return ((cost_t) RATIO)*e->opt_cost - e->alg_cost;
}

cost_t edge_weight_three(adv_outedge *e) {
    return ((cost_t) 3) *e->opt_cost - e->alg_cost;
}

cost_t edge_weight_four(adv_outedge *e) {
    return ((cost_t) 4) *e->opt_cost - e->alg_cost;
}

adv_outedge * locate_edge(adversary_vertex *from, adversary_vertex *to) {
    for (auto e: from->edgelist) {
        if (e->to == to) {
            return e;
        }
    }
    return nullptr;
}

void print_vertex_sequence(std::vector<long int> seq) {
    for (int counter = 0; counter < seq.size(); counter++) {
        fprintf(stderr, "Vertex %d/%zu:\n", counter, seq.size());
        adversary_vertex *v = g.get_vert(seq[counter]);
        v->print(stderr);
        print_memory_info(v->mem);
        if (counter < seq.size() - 1) {
            adversary_vertex *vnext = g.get_vert(seq[counter+1]);
            adv_outedge *e = locate_edge(v, vnext);
            e->print(stderr);
        }
    }
}

void bellman_ford() {
    long unsigned int n = factorial(LISTSIZE)*(max_memory+1);
    fprintf(stderr, "There are %ld vertices in the graph.\n", n);

    cost_t distances[n];
    long int pred[n];

    for (long unsigned int i = 0; i < n; i ++) {
        distances[i] = (cost_t) INT64_MAX;
        pred[i] = -1;
    }

    distances[0] = 0;
    pred[0] = 0;

    for (int iteration = 0; iteration < n; iteration++) {
        fprintf(stderr, "Iteration %d.\n", iteration);
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
                    reverse(cycle.begin(), cycle.end());
                    print_vertex_sequence(cycle);
                    return;
                }
            }
        }
    }
    fprintf(stderr, "No negative cycles present.\n");
    return;
}