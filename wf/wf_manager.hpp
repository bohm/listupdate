#pragma once

#include <random>
#include <unordered_map>
#include "perm_manager.hpp"
#include "workfunction.hpp"


// Mersenne twister.
std::mt19937_64 gen(12345);

uint64_t rand_64bit() {
    uint64_t r = gen();
    return r;
}

template <int SIZE> class wf_manager {
public:
    perm_manager<SIZE>& pm;

    std::array<std::array<uint64_t, diameter_bound(SIZE)+1>, factorial(SIZE)> zobrist;

    std::vector<workfunction<SIZE>> reachable_wfs;
    std::vector<std::array<unsigned int, SIZE>> adjacent_functions;
    std::vector<std::array<short, SIZE>> min_update_costs;
    std::unordered_map<uint64_t, unsigned int> hash_to_index;
    wf_manager(perm_manager<SIZE> &p) : pm(p) {

        for (int i = 0; i < factorial(SIZE); i++) {
            for (int v = 0; v < diameter_bound(SIZE)+1; v++) {
                zobrist[i][v] = rand_64bit();
            }
        }

        initialize_inversions();
    }

    void initialize_inversions() {
        invs.vals[0] = 0;
        for (int i = 1; i < factorial(TESTSIZE); i++) {
            invs.vals[i] = diameter_bound(TESTSIZE);
        }

        dynamic_update(&invs);
    }

    void flat_update(workfunction<SIZE> *wf, short req) {
        for (int i = 0; i < factorial(SIZE); i++) {
            wf->vals[i] += pm.all_perms[i].position(req);
        }
    }

    void cut_minimum(workfunction<SIZE> *wf) {
        short m = wf->min();
        for (int i = 0; i < factorial(SIZE); i++) {
            wf->vals[i] -= m;
        }
    }

    uint64_t hash(workfunction<SIZE> *wf) {
        uint64_t ret = 0;
        for (int i = 0; i < factorial(SIZE); i++) {
            ret ^= zobrist[i][wf->vals[i]];
        }
        return ret;
    }

    void dynamic_update(workfunction<SIZE> *wf) {
        for (short value = 0; value < diameter_bound(SIZE); value++) {
            for (int i = 0; i < factorial(SIZE); i++) {
                if (wf->vals[i] == value) {
                    for (uint64_t adj: pm.adjacencies[i]) {
                        wf->vals[adj] = std::min((short) (value+1), wf->vals[adj]);
                    }
                }
            }
        }
    }

    // Old version.

    uint64_t adjacency_explicit(uint64_t wf_index, short request) {
        workfunction<SIZE> new_wf = reachable_wfs[wf_index];
        flat_update(&new_wf, request);
        cut_minimum(&new_wf);
        dynamic_update(&new_wf);
        // new_wf.validate();
        uint64_t new_wf_index = hash_to_index[hash(&new_wf)];
        return new_wf_index;
    }

    uint64_t adjacency(uint64_t wf_index, short request) {
        return adjacent_functions[wf_index][request];
    }

    short update_cost(uint64_t wf_index, short request) {
        return min_update_costs[wf_index][request];
    }

    void initialize_reachable() {
        std::unordered_set<uint64_t> reachable_hashes;
        std::vector<std::array<uint64_t, SIZE>> adjacencies_by_hash;

        workfunction<SIZE> initial = invs;
        std::queue<workfunction<SIZE>> q;
        reachable_hashes.insert(hash(&initial));
        q.push(initial);
        while (!q.empty()) {
            workfunction<SIZE> front = q.front();
            q.pop();
            hash_to_index[hash(&front)] = reachable_wfs.size();
            reachable_wfs.push_back(front);
            std::array<uint64_t, SIZE> adj;
            std::array<short, SIZE> upd_cost;

            for (short req = 0; req < SIZE; req++) {
                workfunction<SIZE> new_wf = front;
                flat_update(&new_wf, req);
                upd_cost[req] = new_wf.min();
                cut_minimum(&new_wf);
                dynamic_update(&new_wf);
                // new_wf.validate();

                uint64_t h = hash(&new_wf);
                adj[req] = h;
                if (!reachable_hashes.contains(h)) {
                    reachable_hashes.insert(h);
                    q.push(new_wf);
                }
            }
            adjacencies_by_hash.push_back(adj);
            min_update_costs.push_back(upd_cost);
        }

        // Recompute adjacencies to refer to indices, not hashes.
        for (const auto adj_by_hash: adjacencies_by_hash) {
            std::array<unsigned int, SIZE> adj_by_index;
            for (int i = 0; i < SIZE; i++) {
                adj_by_index[i] = hash_to_index[adj_by_hash[i]];
            }
            adjacent_functions.push_back(adj_by_index);
        }

        fprintf(stderr, "Sanity check: set %zu, vector %zu, map %zu, adjacency lists %zu,"
                        " update costs %zu.\n", reachable_hashes.size(),
                reachable_wfs.size(), hash_to_index.size(), adjacent_functions.size(), min_update_costs.size());
    }

};
