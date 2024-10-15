#pragma once

#include "graph_bipartite_mru.hpp"

template <int SIZE> void add_to(std::array<int, SIZE>& one, const std::array<int, SIZE>& two) {
    for (int i = 0; i < SIZE; i++) {
        one[i] += two[i];
    }
}

class bip_mru_potentials {
public:
    graph_bipartite_mru *g;

    static short transform_request(short request, permutation<LISTSIZE> *transformation) {
        return transformation->data[request];
    }

    std::array<std::array<uint64_t, factorial[LISTSIZE]>, factorial[LISTSIZE] * LISTSIZE> *equivalence_classes;

    void compute_equivalence_classes() {
        equivalence_classes = new std::array<std::array<uint64_t, factorial[LISTSIZE]>, factorial[LISTSIZE] * LISTSIZE>;
#pragma omp parallel for
        for (uint64_t mru_index = 0; mru_index < factorial[LISTSIZE]; mru_index++) {
            for (short req = 0; req < LISTSIZE; req++) {
                uint64_t equivalence_index = mru_index * LISTSIZE + req;

                for (uint64_t transformation = 0; transformation < factorial[LISTSIZE]; transformation++) {
                    memory_perm new_mem(mru_index);
                    new_mem = new_mem.recompute(&(g->pm->all_perms[transformation].data));
                    uint64_t equivalent_position = graph_bipartite_mru::alg_index(transformation, new_mem.data,
                                                                                  transform_request(req,
                                                                                                    &(g->pm->all_perms[transformation])));
                    (*equivalence_classes)[equivalence_index][transformation] = equivalent_position;
                }
            }
        }
    }

    void print_equivalence_class(uint64_t eq_class) const {
        for (int i = 0; i < factorial[LISTSIZE]; i++) {
            uint64_t alg_index = (*equivalence_classes)[eq_class][i];
            g->alg_verts[alg_index]->print(stderr);
        }
    }

    bool update_adv() {
        bool any_potential_changed = false;
// #pragma omp parallel for
        for (uint64_t index = 0; index < g->adv_verts.size(); index++) {
            gbm_adversary_vertex *v = g->adv_verts[index];
            if (GRAPH_DEBUG) {
                fprintf(stderr, "ADV vertex update of ");
                v->print(stderr);
            }
            // workfunction<SIZE> wf_before_move = wf.reachable_wfs[wf_index];
            cost_t new_pot = std::numeric_limits<cost_t>::lowest();

            // fprintf(stderr, "Number of edges: %zu.\n", v->edgelist.size());
            for (auto *e: v->edgelist) {
                // fprintf(stderr, "Potential calculation %f vs new_pot %f for edge ", e->to->alg_potential - e->opt_cost, new_pot);
                // e->print(stderr);
                if (e->to->alg_potential - e->opt_cost > new_pot) {
                    if (GRAPH_DEBUG) {
                        fprintf(stderr, "For vertex adv%lu, a potentially improving edge: ", index);
                        e->print(stderr);
                        fprintf(stderr, "Leading to vertex with potential %f: ", e->to->alg_potential);
                        e->to->print(stderr);
                    }
                    new_pot = e->to->alg_potential - e->opt_cost;
                }
            }

            if (v->adv_potential != new_pot) {
                any_potential_changed = true;
                if (GRAPH_DEBUG) {
                    fprintf(stderr, "ADV vertex %" PRIu64 " changed its potential from %f to %f.\n",
                            index, v->adv_potential, new_pot);
                }
                v->adv_potential = new_pot;
            }
        }
        return any_potential_changed;
    }

    bool update_alg() {
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < g->alg_verts.size(); index++) {
            gbm_algorithm_vertex *v = g->alg_verts[index];
            cost_t new_pot = std::numeric_limits<cost_t>::max();

            for (auto *e: v->edgelist) {
                if (e->to->adv_potential + e->alg_cost < new_pot) {
                    new_pot = e->to->adv_potential + e->alg_cost;
                }
            }

            if (v->alg_potential != new_pot) {
                any_potential_changed = true;
                if (GRAPH_DEBUG) {
                    fprintf(stderr, "ALG vertex %" PRIu64 " changed its potential from %f to %f.\n",
                            index, v->alg_potential, new_pot);
                }
                v->alg_potential = new_pot;
            }
        }
        return any_potential_changed;
    }

    bool update_alg_equivalence_classes() const {
        bool any_potential_changed = false;
// #pragma omp parallel for
        for (uint64_t equivalence_class = 0; equivalence_class < factorial[LISTSIZE] * LISTSIZE; equivalence_class++) {
            auto req = (short) (equivalence_class % LISTSIZE);
            uint64_t mru = equivalence_class / LISTSIZE;
            // For the representative of this equivalence class, the request req is at position req.
            // cost_t min_max_potential = std::numeric_limits<cost_t>::max();
            cost_t min_sum_potential = std::numeric_limits<cost_t>::max();
            short minimizer_shift = 0;

            if (GRAPH_DEBUG) {
                fprintf(stderr, "Equivalence class %lu before:\n", equivalence_class);
                print_equivalence_class(equivalence_class);
            }

            // The order of edgelist is the reverse order, so it is the order of "how many places we should shift by.

            for (short shift_by = 0; shift_by <= req; shift_by++) {
                if (GRAPH_DEBUG) {
                    fprintf(stderr, "Considering shift of request by %hd places.\n", shift_by);
                }
                // cost_t max_potential_for_fixed_move = std::numeric_limits<cost_t>::lowest();
                cost_t sum_potential_for_fixed_move = 0;
                for (uint64_t i = 0; i < factorial[LISTSIZE]; i++) {
                    uint64_t alg_index = (*equivalence_classes)[equivalence_class][i];
                    gbm_algorithm_vertex *v = g->get_alg_vert(alg_index);
                    gbm_alg_outedge *e = v->edgelist[shift_by];
                    if (GRAPH_DEBUG) {
                        fprintf(stderr, "For vertex alg%lu, this means edge: ", alg_index);
                        e->print(stderr);
                        fprintf(stderr, "Leading to vertex: with potential %f: ", e->to->adv_potential);
                        e->to->print(stderr);
                    }
                    sum_potential_for_fixed_move += e->to->adv_potential + e->alg_cost;

                    /*
                    if (e->to->adv_potential + e->alg_cost > max_potential_for_fixed_move) {
                        max_potential_for_fixed_move = e->to->adv_potential + e->alg_cost;
                    }
                    */
                }

                if (sum_potential_for_fixed_move < min_sum_potential) {
                    min_sum_potential = sum_potential_for_fixed_move;
                    minimizer_shift = shift_by;
                }

                /*
                if (max_potential_for_fixed_move < min_max_potential) {
                    min_max_potential = max_potential_for_fixed_move;
                    minimizer_shift = shift_by;
                }*/
            }
            fprintf(stderr, "Minimizer shift for class %lu is %hd, with value %f.\n", equivalence_class,
                    minimizer_shift, min_sum_potential);
            // Now, update the potential of all vertices by the minimizer move.
            for (uint64_t i = 0; i < factorial[LISTSIZE]; i++) {
                uint64_t alg_index = (*equivalence_classes)[equivalence_class][i];
                gbm_algorithm_vertex *v = g->get_alg_vert(alg_index);
                gbm_alg_outedge *e = v->edgelist[minimizer_shift];
                cost_t new_pot = e->to->adv_potential + e->alg_cost;
                if (v->alg_potential != new_pot) {
                    any_potential_changed = true;
                    v->alg_potential = new_pot;
                }
            }

            if (GRAPH_DEBUG) {
                fprintf(stderr, "Equivalence class %lu after:\n", equivalence_class);
                print_equivalence_class(equivalence_class);
            }

        }
        return any_potential_changed;
    }


    bool update_alg_frequency_of_min() const {
        bool any_potential_changed = false;
        for (uint64_t equivalence_class = 0; equivalence_class < factorial[LISTSIZE] * LISTSIZE; equivalence_class++) {
            auto req = (short) (equivalence_class % LISTSIZE);
            uint64_t mru = equivalence_class / LISTSIZE;
            std::array<int, LISTSIZE> shift_min_frequency = {};

            if (GRAPH_DEBUG) {
                fprintf(stderr, "Equivalence class %lu before:\n", equivalence_class);
                print_equivalence_class(equivalence_class);
            }

            // The order of edgelist is the reverse order, so it is the order of "how many places we should shift by.
            for (uint64_t i = 0; i < factorial[LISTSIZE]; i++) {
                std::array<int, LISTSIZE> local_min_move = {0};
                cost_t min_pot = std::numeric_limits<cost_t>::max();
                // First phase: Compute the minimum potential.
                for (short shift_by = 0; shift_by <= req; shift_by++) {
                    if (GRAPH_DEBUG) {
                        fprintf(stderr, "Considering shift of request by %hd places.\n", shift_by);
                    }
                    uint64_t alg_index = (*equivalence_classes)[equivalence_class][i];
                    gbm_algorithm_vertex *v = g->get_alg_vert(alg_index);
                    gbm_alg_outedge *e = v->edgelist[shift_by];
                    if (GRAPH_DEBUG) {
                        fprintf(stderr, "For vertex alg%lu, this means edge: ", alg_index);
                        e->print(stderr);
                        fprintf(stderr, "Leading to vertex: with potential %f: ", e->to->adv_potential);
                        e->to->print(stderr);
                    }
                    if (e->to->adv_potential + e->alg_cost < min_pot) {
                        min_pot = e->to->adv_potential + e->alg_cost;
                    }
                }
                // Second phase, mark all shifts with the minimum value as being the minimizers.
                for (short shift_by = 0; shift_by <= req; shift_by++) {
                    uint64_t alg_index = (*equivalence_classes)[equivalence_class][i];
                    gbm_algorithm_vertex *v = g->get_alg_vert(alg_index);
                    gbm_alg_outedge *e = v->edgelist[shift_by];
                    if (e->to->adv_potential + e->alg_cost == min_pot) {
                        local_min_move[shift_by]++;
                    }
                }
                // And add your local min moves to the global frequency array.
                add_to<LISTSIZE>(shift_min_frequency, local_min_move);
            }

            int max_frequency_value = *std::max_element(std::begin(shift_min_frequency), std::end(shift_min_frequency));

            int max_first_position = 0;
            for (;max_first_position < LISTSIZE; max_first_position++) {
                if (shift_min_frequency[max_first_position] == max_frequency_value) {
                    break;
                }
            }

            fprintf(stderr, "The move that is most often the minimum frequency for class %lu is %hd, with value %d.\n",
                    equivalence_class, max_first_position, shift_min_frequency[max_first_position]);
            // Now, update the potential of all vertices by the minimizer move.
            for (uint64_t i = 0; i < factorial[LISTSIZE]; i++) {
                uint64_t alg_index = (*equivalence_classes)[equivalence_class][i];
                gbm_algorithm_vertex *v = g->get_alg_vert(alg_index);
                gbm_alg_outedge *e = v->edgelist[max_first_position];
                cost_t new_pot = e->to->adv_potential + e->alg_cost;
                if (v->alg_potential != new_pot) {
                    any_potential_changed = true;
                    v->alg_potential = new_pot;
                }
            }

            if (GRAPH_DEBUG) {
                fprintf(stderr, "Equivalence class %lu after:\n", equivalence_class);
                print_equivalence_class(equivalence_class);
            }

        }
        return any_potential_changed;
    }
};