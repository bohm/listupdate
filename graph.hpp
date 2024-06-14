#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include "memory_pairs.hpp"
#include "permutations.hpp"
#include "algorithm.hpp"

class adversary_vertex;

using vert_container = std::array<adversary_vertex*, MEMORY::max + 1>;

class graph;

class graph {
public:
    uint64_t edgecounter = 0;
    uint64_t reachable_vertices = 0;
    std::array<vert_container, factorial(LISTSIZE)> verts;

    adversary_vertex *get_vert(permutation *perm, MEMORY m) const
    {
        return verts[lexindex_quadratic(perm)][m.data];
    }

    adversary_vertex *get_vert(long int id);

    void dfs_reachability();
    void reachability_recursive(adversary_vertex *v);
    void reachability_nonrecursive(adversary_vertex *start);
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
    MEMORY mem;
    std::vector<adv_outedge*> edgelist = {};
    bool reachable = false;


    adversary_vertex(permutation *p, MEMORY m) {
        perm = *p;
        mem = m;
        id = lexindex_quadratic(&perm) * (MEMORY::max + 1) + mem.data;
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
            MEMORY mem_copy(mem);
            int alg_cost = ALG_SINGLE_STEP(&perm_copy, &mem_copy, item);
            int opt_cost = item;
            if (FRONT_ACCESS_COSTS_ONE) {
                opt_cost += 1;
            }
            adversary_vertex *target = g.get_vert(&perm_copy, mem_copy);
            auto *edge = new adv_outedge(g.edgecounter++, this, target, item, alg_cost, opt_cost);
            edgelist.push_back(edge);


            // These are purely for testing purposes. We should move them elsewhere ultimately.
            // auto [targid, cost] = implicit_graph::presentation_edge(perm, mem, item);
            // assert(targid == target->id && cost == EDGE_WEIGHT(opt_cost, alg_cost));
        }
    }

    void build_translation_edges() {
        for (int opt_swap = 0; opt_swap < LISTSIZE-1; opt_swap++) {
            permutation single_swap = IDENTITY;
            swap(&single_swap, opt_swap);

            permutation perm_copy(perm);
            MEMORY mem_copy = mem.recompute(&single_swap);
            recompute_alg_perm(&perm_copy, &single_swap);
            adversary_vertex *target = g.get_vert(&perm_copy, mem_copy);
            auto *edge = new adv_outedge(g.edgecounter++, this, target, opt_swap);
            edgelist.push_back(edge);

            // These are purely for testing purposes. We should move them elsewhere ultimately.
            // auto [targid, cost] = implicit_graph::translation_edge(perm, mem, opt_swap);
            // assert(targid == target->id && cost == EDGE_WEIGHT(1, 0));
        }
    }


    void print(FILE *f) {

        fprintf(f, "%lu [label=\"%lu,", id, mem.data);
        print_permutation(&perm, f, false);
        fprintf(f, "\"];\n");

    }

    void print_edges(FILE *f) {
        for (auto e: edgelist) {
            e->print(f);
        }
    }


};


void adv_outedge::print(FILE* f) {
    if (presented_item == -1) {
        fprintf(f, "%lu -> %lu [label=\"swap %d,%d\"];\n", from->id, to->id, opt_swap, opt_swap+1);
    } else {
        fprintf(f, "%lu -> %lu [label=\"req: %d, a_cost: %f, o_cost: %f\"];\n", from->id, to->id, presented_item,
                alg_cost, opt_cost);
    }

}

void add_vertex_to_graph(permutation *perm, MEMORY m) {
    auto *v = new adversary_vertex(perm, m);
    g.verts[lexindex_quadratic(perm)][m.data] = v;
}


adversary_vertex* graph::get_vert(long int id) {
        uint64_t memory_section = id % (MEMORY::max + 1);
        uint64_t permutation_section = id / (MEMORY::max + 1);
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

cost_t edge_weight_param(adv_outedge *e) {
    return edge_weight_param(e->opt_cost, e->alg_cost);
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
        fprintf(stderr, "Memory content for vertex %d/%zu:\n", counter, seq.size());
        v->mem.full_print();
        fprintf(stderr, "\n");
        if (counter < seq.size() - 1) {
            adversary_vertex *vnext = g.get_vert(seq[counter+1]);
            adv_outedge *e = locate_edge(v, vnext);
            e->print(stderr);
        }
    }
}

void graph::dfs_reachability() {
    // reachability_recursive(get_vert(0));
    reachability_nonrecursive(get_vert(0));
    fprintf(stderr, "Graph: %" PRIu64 " vertices were reachable.\n", reachable_vertices);
}

void graph::reachability_nonrecursive(adversary_vertex *start) {
    std::unordered_set<adversary_vertex*> set;
    std::unordered_set<adversary_vertex*> visited;
    set.insert(start);

    while(!set.empty()) {
        adversary_vertex* v = *(set.begin());
        set.erase(v);
        v->reachable = true;
        reachable_vertices++;
        visited.insert(v);
        for (auto& e: v->edgelist) {
            if (!visited.contains(e->to)) {
                set.insert(e->to);
            }
        }
    }
}

void graph::reachability_recursive(adversary_vertex *v) {
    if (v->reachable) {
        return;
    } else {
        v->reachable = true;
        reachable_vertices++;
        for (auto& e: v->edgelist) {
            reachability_recursive(e->to);
        }
    }
}
