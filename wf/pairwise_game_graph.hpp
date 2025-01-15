#pragma once

#pragma once

#include <omp.h>
#include "../pairwise_wf_manager.hpp"
#include "../permutation_graph.hpp"

template <short SIZE> class pairwise_game_graph {
private:
    bool *alg_vertices_visited = nullptr;
    bool *adv_vertices_visited = nullptr;
public:
    short *alg_vertices = nullptr;
    uint64_t algsize = 0;
    short *adv_vertices = nullptr;
    uint64_t advsize = 0;
    pairwise_wf_manager<SIZE> &wf;
    permutation_graph<SIZE> &pg;

    pairwise_game_graph(pairwise_wf_manager<SIZE> &w, permutation_graph<SIZE> &p) : wf(w), pg(p) {
        advsize = wf.reachable_wfs.size() * factorial[SIZE];
        adv_vertices = new short[advsize];
        algsize = wf.reachable_wfs.size() * factorial[SIZE] * SIZE;
        alg_vertices = new short[algsize];

        for (int i = 0; i < algsize; i++) {
            alg_vertices[i] = 0;
        }

        for (int i = 0; i < advsize; i++) {
            adv_vertices[i] = 0;
        }

    }

    ~pairwise_game_graph()
    {
        delete[] adv_vertices;
        delete[] alg_vertices;
    }


    [[nodiscard]] std::pair<unsigned long int, unsigned long int> decode_adv(uint64_t index) const {
        return {index / factorial[SIZE], index % factorial[SIZE]};
    }

    uint64_t encode_adv(unsigned long int wf_index, unsigned long int perm_index) const {
        return wf_index * factorial[SIZE] + perm_index;
    }

    [[nodiscard]] std::tuple<unsigned long int, unsigned long int, unsigned long int> decode_alg(uint64_t index) const {
        unsigned long int request_index = index % SIZE;
        unsigned long int rest = index / SIZE;
        unsigned long int perm_index = rest % factorial[SIZE];
        unsigned long int wf_index = rest / factorial[SIZE];
        return {wf_index, perm_index, request_index};
    }

    uint64_t encode_alg(unsigned long int wf_index, unsigned long int perm_index, unsigned long int request_index) {
        return wf_index * factorial[SIZE] * SIZE + perm_index * SIZE + request_index;
    }

    void print_adv(uint64_t index) const {
        auto [wf_index, perm_index] = decode_adv(index);
        fprintf(stderr, "WF index %lu, perm_index %lu.\n", wf_index, perm_index);
        pg.all_perms[perm_index].print();
        fprintf(stderr, "Work function %lu: \n", wf_index);
        // wf.reachable_wfs[wf_index].print();
    }

    void print_alg(uint64_t index) const {
        auto [wf_index, perm_index, req] = decode_alg(index);
        fprintf(stderr, "ALG vertex: index %lu, perm_index %lu, request %lu.\n", wf_index, perm_index, req);
        pg.all_perms[perm_index].print();
        fprintf(stderr, "Work function %lu: \n", wf_index);
        // wf.reachable_wfs[wf_index].print();
    }

    short adv_cost(unsigned int wf_index, short req) {
        return RATIO*ADV_MULTIPLIER*wf.update_cost(wf_index, req);
    }

    short alg_cost(unsigned int perm_index_one, unsigned int perm_index_two, short req) {
        const permutation<SIZE>& perm_one = pg.all_perms[perm_index_one];
        const permutation<SIZE>* perm_two = &(pg.all_perms[perm_index_two]);
        return ALG_MULTIPLIER*(perm_one.position(req) + perm_one.inversions_wrt(perm_two));
    }

    short min_adv_potential() {
        short m = std::numeric_limits<short>::max();
        for (uint64_t index = 0; index < advsize; index++) {
            m = std::min(m, adv_vertices[index]);
        }
        return m;
    }

    bool update_adv() {
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < advsize; index++) {
            auto [wf_index, perm_index] = decode_adv(index);
            if(GRAPH_DEBUG) {
                fprintf(stderr, "ADV vertex update %" PRIu64 " corresponding to wf_index %lu and perm_index %lu.\n",
                        index, wf_index, perm_index);
                print_adv(index);
            }
            // workfunction<SIZE> wf_before_move = wf.reachable_wfs[wf_index];
            short new_pot = std::numeric_limits<short>::min();
            for (int r = 0; r < SIZE; r++) {
                uint64_t new_wf_index = wf.adjacency(wf_index, r);
                uint64_t alg_index = encode_alg(new_wf_index, perm_index, r);
                short adv_cost_s = adv_cost(wf_index, r);
                if(GRAPH_DEBUG) {
                    fprintf(stderr, "ADV vertex %" PRIu64 " has request %d associated with cost %hd.\n",
                            index, r, adv_cost_s);
                    // ADV potential update.
                    fprintf(stderr, "Phi_x (%hd) - d_yx (%hd) = %hd.\n",
                            alg_vertices[alg_index], adv_cost_s,
                            alg_vertices[alg_index] - adv_cost_s);
                }
                if (alg_vertices[alg_index] - adv_cost_s > new_pot) {
                    new_pot = alg_vertices[alg_index] - adv_cost_s;
                }
            }

            if (adv_vertices[index] != new_pot) {
                any_potential_changed = true;
                if(GRAPH_DEBUG) {
                    fprintf(stderr, "ADV vertex %" PRIu64 " changed its potential from %hd to %hd.\n",
                            index, adv_vertices[index], new_pot);
                }
                adv_vertices[index] = new_pot;
            }
        }
        return any_potential_changed;
    }

    bool update_alg() {
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < algsize; index++) {
            auto [wf_index, perm_index, req] = decode_alg(index);
            if(GRAPH_DEBUG) {
                fprintf(stderr, "ALG vertex update %" PRIu64 " corresponding to wf_index %lu, perm_index %lu, request "
                                "%lu.\n",
                        index, wf_index, perm_index, req);

                print_alg(index);
            }
            short new_pot = std::numeric_limits<short>::max();

            for (int p = 0; p < factorial[SIZE]; p++) {
                uint64_t target_adv = encode_adv(wf_index, p);
                short alg_cost_s = alg_cost(perm_index, p, req);

                // ALG potential update.
                if(GRAPH_DEBUG) {
                    fprintf(stderr, "Phi_y (%hd) + c_xy  (%hd) = %hd.\n",
                            adv_vertices[target_adv], alg_cost_s, adv_vertices[target_adv] + alg_cost_s);
                }
                if (adv_vertices[target_adv] + alg_cost_s < new_pot) {
                    new_pot = adv_vertices[target_adv] + alg_cost_s;
                }
            }

            if (alg_vertices[index] != new_pot) {
                any_potential_changed = true;
                if(GRAPH_DEBUG) {
                    fprintf(stderr, "ALG vertex %" PRIu64 " changed its potential from %hd to %hd.\n",
                            index, alg_vertices[index], new_pot);
                }
                alg_vertices[index] = new_pot;

            }
        }

        return any_potential_changed;
    }


    void print_shortest_path(unsigned long index_adv, const std::vector<short>& shrp) {
        fprintf(stderr, "adv%lu shortest request path: ", index_adv);
        for (auto x: shrp) {
            fprintf(stderr, "%hd, ", x);
        }
        fprintf(stderr, "\n");
    }

};