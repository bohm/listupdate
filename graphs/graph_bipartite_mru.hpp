#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include "../memory_pairs.hpp"
#include "../old_perm_functions.hpp"
#include "../algorithm.hpp"
#include "../wf_manager.hpp"

/* The graph in this case is:
 * bipartite (so n!*r moves for OPT, OPT can swap to any permutation but switches to an alg vertex)
 * stores ALG position (OPT's list always [0 1 2 3 4]) and MRU values.
 * contains multiple choices of moves for ALG.
 */

class gbm_adversary_vertex;

class gbm_algorithm_vertex;

using vert_container = std::array<gbm_adversary_vertex *, memory_perm::max + 1>;

class graph_bipartite_mru;

class graph_bipartite_mru {
public:
    permutation_graph<LISTSIZE> *pm = nullptr;
    uint64_t edgecounter = 0;
    uint64_t reachable_vertices = 0;
    static constexpr uint64_t adv_vertex_range = factorial[LISTSIZE] * factorial[LISTSIZE];
    static constexpr uint64_t alg_vertex_range = factorial[LISTSIZE] * factorial[LISTSIZE] * LISTSIZE;

    static inline uint64_t adv_index(uint64_t alg_position_index, uint64_t mru_index) {
        return alg_position_index * factorial[LISTSIZE] + mru_index;
    }

    static inline std::pair<uint64_t, uint64_t> adv_deindex(uint64_t adv_index) {
        return {adv_index / factorial[LISTSIZE], adv_index % factorial[LISTSIZE]};
    }

    static inline uint64_t alg_index(uint64_t alg_position_index, uint64_t mru_index, short req) {
        return alg_position_index * (factorial[LISTSIZE] * LISTSIZE) + mru_index * LISTSIZE + req;
    }

    static inline std::tuple<uint64_t, uint64_t, short> alg_deindex(uint64_t alg_index) {
        uint64_t alg_position_index_part = alg_index / (factorial[LISTSIZE] * LISTSIZE);
        uint64_t intermediate = alg_index % (factorial[LISTSIZE] * LISTSIZE);
        return {alg_position_index_part, intermediate / LISTSIZE, intermediate % LISTSIZE};
    }

    std::array<gbm_adversary_vertex *, adv_vertex_range> adv_verts;
    std::array<gbm_algorithm_vertex *, alg_vertex_range> alg_verts;


    graph_bipartite_mru(permutation_graph<LISTSIZE> *p) {
        pm = p;
        /*
        if (!inversions_ready) {
            wf_manager<LISTSIZE>::initialize_inversions();
        }
        */
    }

    gbm_adversary_vertex *get_adv_vert(array_as_permutation *perm, memory_perm m) const {
        return adv_verts[adv_index(lexindex_quadratic(perm), m.data)];
    }

    gbm_adversary_vertex *get_adv_vert(uint64_t id) {
        return adv_verts[id];
    }

    gbm_algorithm_vertex *get_alg_vert(uint64_t id) {
        return alg_verts[id];
    }

    void dfs_reachability();

    void reachability_nonrecursive(gbm_adversary_vertex *start);

    void add_adv_vertex(array_as_permutation *perm, memory_perm m);

    void add_alg_vertex(array_as_permutation *perm, memory_perm m, short req);

    void add_adv_outedge(gbm_adversary_vertex *from, gbm_algorithm_vertex *to, short request, cost_t opt_cost);

    void add_alg_outedge(gbm_algorithm_vertex *from, gbm_adversary_vertex *to, cost_t alg_cost);

    void populate_vertices();

    void add_all_outedges(gbm_adversary_vertex *from);
    // In principle, we could add all outedges for the algorithm as well, but that wouldn't make sense for our use.
    void add_all_outedges() {
        for (uint64_t i = 0; i < adv_verts.size(); i++) {
            add_all_outedges(adv_verts[i]);
        }
    }


    void specific_algorithm_outedge(gbm_algorithm_vertex *from);

    void specific_algorithm_outedges() {
        for (uint64_t j = 0; j < alg_verts.size(); j++) {
            specific_algorithm_outedge(alg_verts[j]);
        }
    }


    cost_t min_adv_potential();


    void request_improving_outedges(gbm_algorithm_vertex *from);
    void request_improving_outedges() {
        for (uint64_t j = 0; j < alg_verts.size(); j++) {
            request_improving_outedges(alg_verts[j]);
        }
    }

    void print(FILE *f);
};

graph_bipartite_mru gbm(pg);

class gbm_adv_outedge {
public:
    gbm_adversary_vertex *from = nullptr;
    gbm_algorithm_vertex *to = nullptr;
    short presented_item = -1;
    cost_t opt_cost = 0;

    // Translation edge.
    gbm_adv_outedge(gbm_adversary_vertex *f, gbm_algorithm_vertex *t, short p_item,
                    cost_t o_cost) {
        from = f;
        to = t;
        presented_item = p_item;
        opt_cost = RATIO*o_cost;
    }

    void print(FILE *f) const;

};

class gbm_alg_outedge {
public:
    gbm_algorithm_vertex *from = nullptr;
    gbm_adversary_vertex *to = nullptr;
    cost_t alg_cost = 0;

    gbm_alg_outedge(gbm_algorithm_vertex *f, gbm_adversary_vertex *t, cost_t a_cost) {
        from = f;
        to = t;
        alg_cost = a_cost;
    }

    void print(FILE *f) const;
};


class gbm_algorithm_vertex {
public:
    uint64_t id;
    array_as_permutation perm;
    memory_perm mem;
    short req;
    std::vector<gbm_alg_outedge *> edgelist{};
    bool reachable = false;
    cost_t alg_potential = 0.0;

    gbm_algorithm_vertex(array_as_permutation *p, memory_perm mru, short r, uint64_t alg_index) {
        perm = *p;
        mem = mru;
        req = r;
        id = alg_index;
    }

    void print(FILE *f) {
        fprintf(f, "alg%lu [ALG=", id);
        print_permutation(&perm, f, false);
        fprintf(f, ", MRU=");
        mem.full_print(f, false);
        fprintf(f, ",req=%hd,pot=%f];\n", req, alg_potential);
    }

    void print_edges(FILE *f) {
        for (auto e: edgelist) {
            e->print(f);
        }
    }
};

// A vertex before OPT presents an item.
class gbm_adversary_vertex {
public:
    uint64_t id;
    array_as_permutation perm;
    memory_perm mem;
    std::vector<gbm_adv_outedge *> edgelist = {};
    bool reachable = false;
    cost_t adv_potential = 0.0;


    gbm_adversary_vertex(array_as_permutation *p, memory_perm m, uint64_t adv_index) {
        perm = *p;
        mem = m;
        id = adv_index;
    }

    inline int position(short item) const {
        for (int i = 0; i < LISTSIZE; i++) {
            if (perm[i] == item) {
                return i;
            }
        }
        return 0;
    }


    void print(FILE *f) {
            fprintf(f, "adv%lu [ALG=", id);
            print_permutation(&perm, f, false);
            fprintf(f, ", MRU=");
            mem.full_print(f, false);
            fprintf(f, ",pot=%f];\n", adv_potential);
    }

    void print_edges(FILE *f) {
        for (auto e: edgelist) {
            e->print(f);
        }
    }


};


void gbm_adv_outedge::print(FILE *f = stderr) const {
    fprintf(f, "adv%lu -> alg%lu: request %hd, ", from->id, to->id, to->req);
    fprintf(f, "relabel from algperm ");
    print_permutation(&(from->perm), f, false);
    fprintf(f, " to ");
    print_permutation(&(to->perm), f, false);
    fprintf(f, ", opt_cost %f.\n", opt_cost);
}


void gbm_alg_outedge::print(FILE *f = stderr) const {
    fprintf(f, "alg%lu -> adv%lu: %hd, ", from->id, to->id, from->req);
    fprintf(f, "update alg from ");
    print_permutation(&(from->perm), f, false);
    fprintf(f, " to ");
    print_permutation(&(to->perm), f, false);
    fprintf(f, ", alg_cost %f.\n", alg_cost);
}

void graph_bipartite_mru::add_adv_vertex(array_as_permutation *perm, MEMORY m) {
    uint64_t adv_index = graph_bipartite_mru::adv_index(lexindex_quadratic(perm), m.data);
    auto *v = new gbm_adversary_vertex(perm, m, adv_index);
    adv_verts[adv_index] = v;
}


void graph_bipartite_mru::add_alg_vertex(array_as_permutation *perm, memory_perm m, short req) {
    uint64_t alg_index = graph_bipartite_mru::alg_index(lexindex_quadratic(perm), m.data, req);
    auto *v = new gbm_algorithm_vertex(perm, m, req, alg_index);
    alg_verts[alg_index] = v;
}


void graph_bipartite_mru::populate_vertices() {
    for (int i = 0; i < factorial[LISTSIZE]; i++) {
        for (int j = 0; j < factorial[LISTSIZE]; j++) {
            permutation<LISTSIZE> &alg_list = pg->all_perms[i];
            memory_perm mru;
            mru.data = j;
            add_adv_vertex(&(pg->all_perms[i].data), mru);


            for (short req = 0; req < LISTSIZE; req++) {
                add_alg_vertex(&(pg->all_perms[i].data), mru, req);
            }
        }
    }
}

void graph_bipartite_mru::print(FILE *f) {
    for (int i = 0; i < adv_verts.size(); i++) {
        adv_verts[i]->print(f);
        adv_verts[i]->print_edges(f);
    }

    for (int j = 0; j < alg_verts.size(); j++) {
        alg_verts[j]->print(f);
        alg_verts[j]->print_edges(f);
    }
}

/*
cost_t edge_weight_param(gbm_adv_outedge *e) {
    return edge_weight_param(e->opt_cost, e->alg_cost);
}
*/

gbm_adv_outedge *locate_advedge(gbm_adversary_vertex *from, gbm_algorithm_vertex *to) {
    for (auto e: from->edgelist) {
        if (e->to == to) {
            return e;
        }
    }
    return nullptr;
}

gbm_alg_outedge *locate_algedge(gbm_algorithm_vertex *from, gbm_adversary_vertex *to) {
    for (auto e: from->edgelist) {
        if (e->to == to) {
            return e;
        }
    }
    return nullptr;
}

void print_vertex_sequence(std::vector<std::pair<bool, long int>> seq) {
    bool is_adv_vertex = false;
    long id = 0;

    for (int counter = 0; counter < seq.size(); counter++) {
        fprintf(stderr, "Vertex %d/%zu:\n", counter, seq.size());
        std::tie(is_adv_vertex, id) = seq[counter];
        if (is_adv_vertex) {
            gbm_adversary_vertex *v = gbm.get_adv_vert(id);
            v->print(stderr);
            fprintf(stderr, "Memory content for vertex %d/%zu:\n", counter, seq.size());
            v->mem.full_print();

            fprintf(stderr, "\n");
            if (counter < seq.size() - 1) {
                std::tie(is_adv_vertex, id) = seq[counter+1];
                gbm_algorithm_vertex *vnext = gbm.get_alg_vert(id);
                gbm_adv_outedge *e = locate_advedge(v, vnext);
                e->print(stderr);
            }
        } else {
            gbm_algorithm_vertex *v = gbm.get_alg_vert(id);
            v->print(stderr);
            fprintf(stderr, "Memory content for vertex %d/%zu:\n", counter, seq.size());
            v->mem.full_print();

            fprintf(stderr, "\n");
            if (counter < seq.size() - 1) {
                std::tie(is_adv_vertex, id) = seq[counter+1];
                gbm_adversary_vertex *vnext = gbm.get_adv_vert(id);
                gbm_alg_outedge *e = locate_algedge(v, vnext);
                e->print(stderr);
            }
        }
    }
}

double total_alg_cost(std::vector<std::pair<bool, long>>  &seq) {
    bool is_adv_vertex = false;
    long id = 0;
    double ret = 0;
    for (int counter = 0; counter < seq.size(); counter++) {
        std::tie(is_adv_vertex, id) = seq[counter];
        if (!is_adv_vertex) {
            gbm_algorithm_vertex *v = gbm.get_alg_vert(id);
            if (counter < seq.size() - 1) {
                std::tie(is_adv_vertex, id) = seq[counter + 1];
                gbm_adversary_vertex *vnext = gbm.get_adv_vert(id);
                gbm_alg_outedge *e = locate_algedge(v, vnext);
                ret += e->alg_cost;
            }
        }
    }
    return ret;
}


double total_opt_cost(std::vector<std::pair<bool, long>> &seq) {
    bool is_adv_vertex = false;
    long id = 0;
    double ret = 0;
    for (int counter = 0; counter < seq.size(); counter++) {
        std::tie(is_adv_vertex, id) = seq[counter];
        if (is_adv_vertex) {
            gbm_adversary_vertex* v = gbm.get_adv_vert(id);
            if (counter < seq.size() - 1) {
                std::tie(is_adv_vertex, id) = seq[counter+1];
                gbm_algorithm_vertex *vnext = gbm.get_alg_vert(id);
                gbm_adv_outedge *e = locate_advedge(v, vnext);
                ret += e->opt_cost;
            }
        }
    }
    return ret;

}

void graph_bipartite_mru::dfs_reachability() {
    // reachability_recursive(get_vert(0));
    reachability_nonrecursive(get_adv_vert(0));
    fprintf(stderr, "Graph: %" PRIu64 " vertices were reachable.\n", reachable_vertices);
}

void graph_bipartite_mru::reachability_nonrecursive(gbm_adversary_vertex *start) {
    std::queue<std::pair<bool, uint64_t>> q;
    std::unordered_set<uint64_t> visited_adv;
    std::unordered_set<uint64_t> visited_alg;

    q.emplace(true, start->id);
    bool cur_adv_v = false;
    uint64_t cur_id = 0;
    while (!q.empty()) {
        std::tie(cur_adv_v, cur_id) = q.front();
        q.pop();
        if (cur_adv_v) {
            adv_verts[cur_id]->reachable = true;
            reachable_vertices++;
            visited_adv.insert(cur_id);
            for (auto &e: adv_verts[cur_id]->edgelist) {
                if (!visited_alg.contains(e->to->id)) {
                    q.emplace(false, e->to->id);
                }
            }
        } else {
            alg_verts[cur_id]->reachable = true;
            reachable_vertices++;
            visited_alg.insert(cur_id);
            for (auto &e: alg_verts[cur_id]->edgelist) {
                if (!visited_adv.contains(e->to->id)) {
                    q.emplace(true, e->to->id);
                }
            }
        }
    }
}

void graph_bipartite_mru::add_adv_outedge(gbm_adversary_vertex *from, gbm_algorithm_vertex *to, short request,
                                          cost_t opt_cost) {
    auto *edge = new gbm_adv_outedge(from, to, request, opt_cost);
    from->edgelist.push_back(edge);
}

void graph_bipartite_mru::add_alg_outedge(gbm_algorithm_vertex *from, gbm_adversary_vertex *to, cost_t alg_cost) {
    auto *edge = new gbm_alg_outedge(from, to, alg_cost);
    from->edgelist.push_back(edge);
}

void graph_bipartite_mru::add_all_outedges(gbm_adversary_vertex *from) {

    for (short newreq = 0; newreq < LISTSIZE; newreq++) {
        for (uint64_t i = 0; i < factorial[LISTSIZE]; i++) {
            permutation<LISTSIZE> opt_relabeling = pm->all_perms[i];
            array_as_permutation perm_after_relabel(from->perm);
            memory_perm mem_after_relabel = from->mem.recompute(&(opt_relabeling.data));
            recompute_alg_perm(&perm_after_relabel, &(opt_relabeling.data));
            uint64_t target_alg_id = alg_index(lexindex_quadratic(&perm_after_relabel), mem_after_relabel.data, newreq);
            gbm_algorithm_vertex *target = get_alg_vert(target_alg_id);
            add_adv_outedge(from, target, newreq, (cost_t) (newreq + opt_relabeling.inversions()));

        }
    }
}

void graph_bipartite_mru::specific_algorithm_outedge(gbm_algorithm_vertex *from) {
    short item = from->req;

    array_as_permutation perm_copy(from->perm);
    memory_perm mem_copy(from->mem);
    cost_t alg_cost = ALG_SINGLE_STEP(&perm_copy, &mem_copy, item);

    gbm_adversary_vertex *target = gbm.get_adv_vert(&perm_copy, mem_copy);

    add_alg_outedge(from, target, alg_cost);
}

cost_t graph_bipartite_mru::min_adv_potential() {
    cost_t m = std::numeric_limits<cost_t>::max();
    for (uint64_t index = 0; index < adv_verts.size(); index++) {
        m = std::min(m, adv_verts[index]->adv_potential);
    }
    return m;
}

void graph_bipartite_mru::request_improving_outedges(gbm_algorithm_vertex *from) {
    short request = from->req;
    memory_perm new_mem = from->mem;
    new_mem.mtf(request);
    permutation<LISTSIZE> perm(from->perm);
    short req_old_position = perm.position(request);
    for (short improving_position = req_old_position; improving_position >= 0; improving_position--) {
        permutation<LISTSIZE> new_perm = perm.move_forward_copy(request, improving_position);
        gbm_adversary_vertex *target = get_adv_vert(&new_perm.data, new_mem);
        cost_t alg_cost = req_old_position + (req_old_position - improving_position);
        add_alg_outedge(from, target, alg_cost);
    }
}

