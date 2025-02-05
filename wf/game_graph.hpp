#pragma once

#include <omp.h>
#include "../wf_manager.hpp"

template <short SIZE> class game_graph {
private:
    bool *alg_vertices_visited = nullptr;
    bool *adv_vertices_visited = nullptr;
public:
    short *alg_vertices = nullptr;
    uint64_t algsize = 0;
    short *adv_vertices = nullptr;
    uint64_t advsize = 0;
    wf_manager<SIZE> &wf;

    bool wfa_adjacencies = false;

    std::vector<uint64_t>* wfa_minimum_adjacency = nullptr;

    game_graph(wf_manager<SIZE> &w, bool wfa_adj = false) : wf(w), wfa_adjacencies(wfa_adj) {
        advsize = wf.reachable_wfs.size()* factorial[SIZE];
        adv_vertices = new short[advsize];
        algsize = wf.reachable_wfs.size() * factorial[SIZE] * SIZE;
        alg_vertices = new short[algsize];

        if (wfa_adjacencies)
        {
            wfa_minimum_adjacency = new std::vector<uint64_t>[algsize];
        }

        for (int i = 0; i < algsize; i++) {
            alg_vertices[i] = 0;
        }

        for (int i = 0; i < advsize; i++) {
            adv_vertices[i] = 0;
        }

    };

    ~game_graph()
    {
        delete[] adv_vertices;
        delete[] alg_vertices;
        if (wfa_adjacencies)
        {
            delete[] wfa_minimum_adjacency;
        }
    }


    std::pair<unsigned long int, unsigned long int> decode_adv(uint64_t index) const {
        return {index / factorial[SIZE], index % factorial[SIZE]};
    }

    uint64_t encode_adv(unsigned long int wf_index, unsigned long int perm_index) const {
        return wf_index * factorial[SIZE] + perm_index;
    }

    std::tuple<unsigned long int, unsigned long int, unsigned long int> decode_alg(uint64_t index) const {
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
        wf.pm.all_perms[perm_index].print();
        fprintf(stderr, "Work function %lu: \n", wf_index);
        // wf.reachable_wfs[wf_index].print();
    }

    void print_alg(uint64_t index) const {
        auto [wf_index, perm_index, req] = decode_alg(index);
        fprintf(stderr, "ALG vertex: index %lu, perm_index %lu, request %lu.\n", wf_index, perm_index, req);
        wf.pm.all_perms[perm_index].print();
        fprintf(stderr, "Work function %lu: \n", wf_index);
        // wf.reachable_wfs[wf_index].print();
    }

    short adv_cost(unsigned int wf_index, short req) {
        return RATIO*MULTIPLIER*wf.update_cost(wf_index, req);
    }

    short alg_cost(unsigned int perm_index_one, unsigned int perm_index_two, short req) {
        const permutation<SIZE>& perm_one = wf.pm.all_perms[perm_index_one];
        const permutation<SIZE>* perm_two = &(wf.pm.all_perms[perm_index_two]);
        return MULTIPLIER*(perm_one.position(req) + perm_one.inversions_wrt(perm_two));
    }

    /*
    short min_adv_potential() {
        short m = std::numeric_limits<short>::max();
        for (uint64_t index = 0; index < advsize; index++) {
            m = std::min(m, adv_vertices[index]);
        }
        return m;
    }
     */

    short min_adv_potential() {
        short m = std::numeric_limits<short>::max();
        for (uint64_t index = 0; index < advsize; index++) {
            m = std::min(m, adv_vertices[index]);
        }
        return m;
    }

    unsigned int wfa_cost(unsigned long wf_index, unsigned long current_alg_index, unsigned long perm_index) const {
        permutation<LISTSIZE>* perm = &(wf.pm.all_perms[perm_index]);
        // permutation<LISTSIZE>* current_alg_pos = &(wf.pm.all_perms[current_alg_index]);
        unsigned int wf_cost = wf.reachable_wfs[wf_index].vals[perm_index];
        // unsigned int transition_cost =  perm->inversions_wrt(current_alg_pos);
        unsigned int quick_transition_cost = wf.pm.quick_inversion_wrt(perm_index, current_alg_index);
        // assert(quick_transition_cost == transition_cost);
        return wf_cost + quick_transition_cost;
    }

    unsigned int workfunction_algorithm_minimum(unsigned long wf_index, unsigned long perm_index, short request) const {
        unsigned int minimum_wfa_cost = std::numeric_limits<unsigned int>::max();
        permutation<LISTSIZE>* current_alg_pos = &(wf.pm.all_perms[perm_index]);
        for (unsigned int i = 0; i < factorial[LISTSIZE]; i++) {
            unsigned int cost_for_i = wfa_cost(wf_index, perm_index, i);
            if (cost_for_i < minimum_wfa_cost) {
                minimum_wfa_cost = cost_for_i;
            }
        }

        return minimum_wfa_cost;
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


    bool update_alg_wfa() {
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < algsize; index++) {
            auto [wf_index, perm_index, req] = decode_alg(index);
            short new_pot = std::numeric_limits<short>::max();

            // Instead of any permutation, we filter those which have higher than minimum value of WFA.
            unsigned int wfa_minimum_value = workfunction_algorithm_minimum(wf_index, perm_index, req);
            for (int p = 0; p < factorial[SIZE]; p++) {
                unsigned int wfa_cost_for_this_index = wfa_cost(wf_index, perm_index, p);
                if (wfa_cost_for_this_index > wfa_minimum_value) {
                    continue;
                }
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


    void build_wfa_adjacencies()
    {
        fprintf(stderr, "Building work function minima adjacencies.\n");
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < algsize; index++)
        {
            auto [wf_index, perm_index, req] = decode_alg(index);
            short new_pot = std::numeric_limits<short>::max();

            // Instead of any permutation, we filter those which have higher than minimum value of WFA.
            unsigned int wfa_minimum_value = workfunction_algorithm_minimum(wf_index, perm_index, req);
            for (int p = 0; p < factorial[SIZE]; p++)
            {
                unsigned int wfa_cost_for_this_index = wfa_cost(wf_index, perm_index, p);
                if (wfa_cost_for_this_index > wfa_minimum_value) {
                    continue;
                } else {
                    uint64_t target_adv = encode_adv(wf_index, p);
                    wfa_minimum_adjacency[index].push_back(target_adv);
                }
            }
        }
        fprintf(stderr, "Adjacencies finalized.\n");
    }

    bool update_alg_wfa_fast() {
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < algsize; index++) {
            auto [wf_index, perm_index, req] = decode_alg(index);
            short new_pot = std::numeric_limits<short>::max();
            for (int p = 0; p < wfa_minimum_adjacency[index].size(); p++) {
                uint64_t target_adv = wfa_minimum_adjacency[index][p];
                auto [wf_index, adv_perm_index] = decode_adv(target_adv);
                short alg_cost_s = alg_cost(perm_index, adv_perm_index, req);
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

    bool update_alg_mtf_of_request() {
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < algsize; index++) {
            auto [wf_index, perm_index, req] = decode_alg(index);
            if(GRAPH_DEBUG) {
                fprintf(stderr, "ALG vertex update %" PRIu64 " corresponding to wf_index %lu, perm_index %lu, request "
                                "%hd.\n",
                        index, wf_index, perm_index, req);

                print_alg(index);
            }
            short new_pot = std::numeric_limits<short>::max();

            // Instead of all choices, we only allow two MTF choices.
            // STAY
            unsigned long stay_perm_index = perm_index;
            uint64_t target_adv = encode_adv(wf_index, stay_perm_index);
            short alg_cost_s = alg_cost(perm_index, stay_perm_index, req);
            if (adv_vertices[target_adv] + alg_cost_s < new_pot) {
                new_pot = adv_vertices[target_adv] + alg_cost_s;
            }

            // MTF

            unsigned long mtf_perm_index = wf.pm.all_perms[perm_index].mtf_copy(req).id();
            target_adv = encode_adv(wf_index, mtf_perm_index);
            alg_cost_s = alg_cost(perm_index, mtf_perm_index, req);
            if (adv_vertices[target_adv] + alg_cost_s < new_pot) {
                new_pot = adv_vertices[target_adv] + alg_cost_s;
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

    bool update_alg_request_moves_forward() {
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

            // Instead of all choices, we only allow the request to move forward.
            short request_pos = wf.pm.all_perms[perm_index].position(req);
            for (short target = 0; target <= request_pos; target++) {
                unsigned long p = wf.pm.all_perms[perm_index].move_forward_copy(req, target).id();


                if(GRAPH_DEBUG) {
                    fprintf(stderr, "alg%lu: Evaluating edge: ", index);
                    wf.pm.all_perms[perm_index].print(stderr, false);
                    fprintf(stderr, " -> ");
                    wf.pm.all_perms[p].print();
                }
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

    bool update_alg_single_swap() {
        bool any_potential_changed = false;
#pragma omp parallel for
        for (uint64_t index = 0; index < algsize; index++) {
            auto [wf_index, perm_index, req] = decode_alg(index);
            if(GRAPH_DEBUG) {
                fprintf(stderr, "ALG vertex update %" PRIu64 " corresponding to wf_index %lu, perm_index %lu, request "
                                "%hd.\n",
                        index, wf_index, perm_index, req);

                print_alg(index);
            }
            short new_pot = std::numeric_limits<short>::max();

            // Instead of all choices, we only allow single swaps.
            for (short swap= 0; swap < SIZE; swap++) {
                unsigned long p = perm_index;
                if (swap < SIZE-1) {
                    p = wf.pm.all_perms[perm_index].swap(swap).id();
                }

                if(GRAPH_DEBUG) {
                    fprintf(stderr, "alg%lu: Evaluating edge: ", index);
                    wf.pm.all_perms[perm_index].print(stderr, false);
                    fprintf(stderr, " -> ");
                    wf.pm.all_perms[p].print();
                }
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

    /*
    short conjectured_potential(unsigned long wf_index, unsigned long perm_index) {
        workfunction<SIZE>& workf = wf.reachable_wfs[wf_index];
        permutation<SIZE> perm = wf.pm.all_perms[perm_index];
        short m = std::numeric_limits<short>::min();
        for (unsigned long p = 0; p < factorial[SIZE]; p++) {
            short inv = perm.inversions_wrt(&(wf.pm.all_perms[p]));
            if (2*inv - 3*workf.vals[p] > m)
            {
                m = 2*inv - 3*workf.vals[p];
            }
        }
        return m;
    }
     */

    void print_shortest_path(unsigned long index_adv, const std::vector<short>& shrp) {
        fprintf(stderr, "adv%lu shortest request path: ", index_adv);
        for (auto x: shrp) {
            fprintf(stderr, "%hd, ", x);
        }
        fprintf(stderr, "\n");
    }



    void print_potential() {
        adv_vertices_visited = new bool[advsize];
        for (int i = 0; i < advsize; i++) {
            adv_vertices_visited[i] = false;
        }
        alg_vertices_visited = new bool[algsize];
        for (int j = 0; j < algsize; j++) {
            alg_vertices_visited[j] = false;
        }

        std::queue<std::tuple<bool, unsigned long, std::vector<short>>> bfs_q;
        std::tuple<bool, unsigned long, std::vector<short>> initial{true, 0, std::vector<short>()};
        bfs_q.emplace(initial);
        adv_vertices_visited[0] = true;
        while(!bfs_q.empty()) {
            auto [adv_vertex, index, shortest_path] = bfs_q.front();
            bfs_q.pop();
            if (adv_vertex) {
                auto [wf_index, perm_index] = decode_adv(index);
                print_adv(index);
                // short conp = conjectured_potential(wf_index, perm_index);
                // fprintf(stderr, "adv%lu: adv potential %hd, marek's conjecture %hd, difference %hd.\n",
                //        index, adv_vertices[index], conp, conp-adv_vertices[index]);
                fprintf(stderr, "adv%lu: adv potential %hd.\n", index, adv_vertices[index]);
                print_shortest_path(index, shortest_path);
                for (short r = 0; r < SIZE; r++) {
                    unsigned long next_wf = wf.adjacency(wf_index, r);
                    unsigned long next_alg_index = encode_alg(next_wf, perm_index, r);
                    std::vector<short> next_shortest_path(shortest_path);
                    next_shortest_path.push_back(r);
                    fprintf(stderr, "adv%lu with req %hd: updated work function number %lu. Next alg%lu.\n",
                            index, r, next_wf, next_alg_index);
                    if (!alg_vertices_visited[next_alg_index]) {
                        alg_vertices_visited[next_alg_index] = true;
                        bfs_q.emplace(false, next_alg_index, next_shortest_path);
                    }
                }
            } else {
                auto [wf_index, perm_index, request] = decode_alg(index);
                print_alg(index);
                fprintf(stderr, "alg%lu: alg potential %hd.\n", index, alg_vertices[index]);

                // Main difference: we only expand on tight edges here.
                // Compute the number of tight edges first. It slows down the printing, but nothing important.
                unsigned long tight_size = 0;
                for (unsigned long p = 0; p < factorial[SIZE]; p++) {
                    unsigned long next_adv_index = encode_adv(wf_index, p);
                    short alg_cost_s = alg_cost(perm_index, p, request);
                    if (adv_vertices[next_adv_index] + alg_cost_s == alg_vertices[index]) {
                        tight_size++;
                    }
                }

                unsigned long tight_index = 0;
                permutation<LISTSIZE> current_alg_pos = wf.pm.all_perms[perm_index];
                unsigned long wfa_minimum = workfunction_algorithm_minimum(wf_index, perm_index, request);
                for (unsigned long p = 0; p < factorial[SIZE]; p++) {
                    unsigned long next_adv_index = encode_adv(wf_index, p);
                    short alg_cost_s = alg_cost(perm_index, p, request);
                    if (adv_vertices[next_adv_index] + alg_cost_s == alg_vertices[index]) {
                        fprintf(stderr, "Tight edge %lu/%lu between alg%lu and adv%lu.",
                                tight_index, tight_size, index, next_adv_index);
                        unsigned long wfa_cost_for_this_move = wfa_cost(wf_index, &current_alg_pos, p);
                        if (wfa_cost_for_this_move == wfa_minimum) {
                            fprintf(stderr, " Minimum of WFA (%lu = %lu).",
                                    wfa_cost_for_this_move, wfa_minimum);
                        } else {
                            fprintf(stderr, " Nonoptimal in terms of WFA (%lu != %lu).",
                                    wfa_cost_for_this_move, wfa_minimum);
                        }

                        if (p == perm_index) {
                            fprintf(stderr, " Stay: ");
                            wf.pm.all_perms[perm_index].print();
                        } else {
                            fprintf(stderr, " Move");
                            if (wf.pm.all_perms[p].data[0] == request) {
                                fprintf(stderr, " (RIF)");
                            }
                            fprintf(stderr, ": ");
                            wf.pm.all_perms[perm_index].print(stderr, false);
                            fprintf(stderr, " -> ");
                            wf.pm.all_perms[p].print();
                        }

                        if (!adv_vertices_visited[next_adv_index]) {
                            adv_vertices_visited[next_adv_index] = true;
                            bfs_q.emplace(true, next_adv_index, shortest_path);
                        }

                        tight_index++;
                    }
                }
            }
        }
        delete[] adv_vertices_visited;
        delete[] alg_vertices_visited;
    }
};