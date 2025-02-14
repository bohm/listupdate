#pragma once

#include <array>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <queue>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <boost/serialization/binary_object.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
using boost::serialization::make_binary_object;


#include "common.hpp"
#include "wf/double_zobrist.hpp"


constexpr uint64_t UNSORTED_PAIRS = (TSIZE * (TSIZE-1))/2;

// Nomenclature:
// value 0: A < B
// value 1: A ?= B
// value 2: A > B
// We also count costs in essentially half-moves, so a cost of 1 is a change from A < B to A ?= B.

template <int SIZE> class pairwise_workfunction {
public:
    std::array<short, UNSORTED_PAIRS> vals;


    // Access a sorted pair, that is, i < j.
    inline int sorted_pair_index(int i, int j) const {
        return canonical_order[i*LISTSIZE+j];
    }

    // Access an unsorted pair, that is, either i > j or i < j (but not i == j).
    inline int unsorted_pair_index(int i, int j) const {
        return sorted_pair_index(std::min(i,j), std::max(i,j));
    }

    void print(FILE *outf = stderr) const {
        for (int i = 0; i < SIZE; i++) {
            for (int j = (i+1); j < SIZE; j++) {
                fprintf(outf, "(%d, %d)->%hd, ", i, j, vals[canonical_order[i * LISTSIZE + j]]);
            }
        }

        fprintf(outf, "\n");
    }

    int update(short request) {
        int cost = 0;
        for (int i = 0; i < request; i++) {
            int index = sorted_pair_index(i, request);
            // Accessing the second element, so we pay 0 only in the case that i > request.
            if (vals[index] <= 1) {
                vals[index]++;
                cost++;
            }
        }

        for (int i = request+1; i < SIZE; i++) {
            int index = sorted_pair_index(request, i);
            // Accessing the first element of the pair, so we pay 0 if the order is 0.
            if (vals[index] > 0) {
                vals[index]--;
                cost++;
            }
        }

        return cost;
    }

    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & vals;
    }
};

template <int SIZE> class pairwise_wf_manager {
public:

    static constexpr pairwise_workfunction<SIZE> initial_partial_workfunction() {
        pairwise_workfunction<SIZE> x;
        for (int i = 0; i < UNSORTED_PAIRS; i++) {
            x.vals[i] = 0;
        }
        return x;
    }
    static constexpr pairwise_workfunction<SIZE> initial_pwf = initial_partial_workfunction();

    std::array<std::array<uint64_t, 3>, UNSORTED_PAIRS> *zobrist;
    std::vector<pairwise_workfunction<SIZE>> reachable_wfs;
    std::unordered_map<uint64_t, unsigned int> hash_to_index;
    std::vector<std::array<unsigned int, SIZE>> adjacent_functions;
    std::vector<std::array<short, SIZE>> update_costs;
    pairwise_wf_manager() {

            zobrist = new std::array<std::array<uint64_t, 3>, UNSORTED_PAIRS>();
            for (int i = 0; i < UNSORTED_PAIRS; i++) {
                for (int v = 0; v < 2; v++) {
                    (*zobrist)[i][v] = rand_64bit();
                }
            }
    }

    ~pairwise_wf_manager() {
        delete zobrist;
    }

    uint64_t hash(pairwise_workfunction<SIZE> *pwf) {
        uint64_t ret = 0;
        for (int i = 0; i < UNSORTED_PAIRS; i++) {
            ret ^= (*zobrist)[i][pwf->vals[i]];
        }
        return ret;
    }


    uint64_t adjacency(uint64_t wf_index, short request) {
        return adjacent_functions[wf_index][request];
    }

    short update_cost(uint64_t wf_index, short request) {
        return update_costs[wf_index][request];
    }


    void save_reachable(std::string reachable_filename) {
        std::ofstream f(reachable_filename.c_str(), std::ofstream::binary);
        boost::archive::binary_oarchive ar(f);

        ar << reachable_wfs;
        ar << adjacent_functions;
        ar << update_costs;
        // ar << make_binary_object(&reachable_wfs, sizeof(reachable_wfs));
        // ar << make_binary_object(hash_to_index, sizeof(hash_to_index));
        // ar << make_binary_object(&adjacent_functions, sizeof(adjacent_functions));
        // ar << make_binary_object(&min_update_costs, sizeof(min_update_costs));
    }

    void load_reachable(std::string reachable_filename) {
        std::ifstream f(reachable_filename.c_str(), std::ifstream::binary);
        boost::archive::binary_iarchive ar(f);
        ar >> reachable_wfs;
        ar >> adjacent_functions;
        ar >> update_costs;
    }

    void fill_hash_to_index() {
        hash_to_index.clear();
        for (int i = 0; i < reachable_wfs.size(); i++) {
            hash_to_index[hash(&(reachable_wfs[i]))] = i;
        }
    }
    void initialize_reachable(std::string reachable_filename) {
        if (std::filesystem::exists(reachable_filename)) {
            load_reachable(reachable_filename);
            fill_hash_to_index();
            fprintf(stderr, "Sanity check: reachable_wfs %zu, hash_to_index %zu, adjacent functions %zu,"
                " update costs %zu.\n", reachable_wfs.size(), hash_to_index.size(), adjacent_functions.size(),
                update_costs.size());

        } else {
            initialize_reachable_from_scratch();
            save_reachable(reachable_filename);
        }
    }


    void initialize_reachable_from_scratch() {
        std::unordered_set<uint64_t> reachable_hashes;

        std::vector<std::array<uint64_t, SIZE>> adjacencies_by_hash;

        pairwise_workfunction<SIZE> initial = initial_pwf;
        std::queue<pairwise_workfunction<SIZE>> q;
        reachable_hashes.insert(hash(&initial));
        q.push(initial);
        while (!q.empty()) {
            pairwise_workfunction<SIZE> front = q.front();
            q.pop();
            hash_to_index[hash(&front)] = reachable_wfs.size();
            reachable_wfs.push_back(front);
            std::array<uint64_t, SIZE> adj;
            std::array<short, SIZE> upd_cost;

            for (short req = 0; req < SIZE; req++) {
                pairwise_workfunction<SIZE> new_wf = front;
                int cost = new_wf.update(req);
                upd_cost[req] = cost;
                // new_wf.validate();

                uint64_t h = hash(&new_wf);
                adj[req] = h;
                if (!reachable_hashes.contains(h)) {
                    reachable_hashes.insert(h);
                    q.push(new_wf);
                }
            }
            adjacencies_by_hash.push_back(adj);
            update_costs.push_back(upd_cost);
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
                reachable_wfs.size(), hash_to_index.size(), adjacent_functions.size(), update_costs.size());
    }


    void print_reachable(const std::string& filename) {
        FILE* outf = fopen(filename.c_str(), "w");
        for (int i = 0;i < reachable_wfs.size(); i++) {
            fprintf(outf, "Reachable pairwise work function of id (index in array) %d:\n", i);
            reachable_wfs[i].print(outf);
        }
        fclose(outf);
    }

};
