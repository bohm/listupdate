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

template <int SIZE> class algorithm_against_pairwise_opt {
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

    static constexpr algorithm_against_pairwise_opt<SIZE> initial_partial_workfunction() {
        algorithm_against_pairwise_opt<SIZE> x;
        for (int i = 0; i < UNSORTED_PAIRS; i++) {
            x.vals[i] = 0;
        }
        return x;
    }
    static constexpr algorithm_against_pairwise_opt<SIZE> initial_pwf = initial_partial_workfunction();

    std::array<std::array<uint64_t, 3>, UNSORTED_PAIRS> *zobrist;
    // std::vector<pairwise_workfunction<SIZE>> reachable_wfs;
    // std::vector<std::array<unsigned int, SIZE>> adjacent_functions;
    // std::vector<std::array<short, SIZE>> update_costs;

    uint64_t reachable_workfunctions = 0;
    algorithm_against_pairwise_opt<SIZE>* reachable_wfs_arr = nullptr;
    std::array<unsigned int, SIZE>* adjacent_functions_arr = nullptr;
    std::array<short, SIZE>* update_costs_arr = nullptr;

    std::unordered_map<uint64_t, unsigned int> hash_to_index;
    // std::vector<std::array<unsigned int, SIZE>> adjacent_functions;
    // std::vector<std::array<short, SIZE>> update_costs;
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

    uint64_t hash(algorithm_against_pairwise_opt<SIZE> *pwf) {
        uint64_t ret = 0;
        for (int i = 0; i < UNSORTED_PAIRS; i++) {
            ret ^= (*zobrist)[i][pwf->vals[i]];
        }
        return ret;
    }


    uint64_t adjacency(uint64_t wf_index, short request) {
        return adjacent_functions_arr[wf_index][request];
    }

    short update_cost(uint64_t wf_index, short request) {
        return update_costs_arr[wf_index][request];
    }

    void serialize_reachable(const std::string& reachable_filename) const {
        FILE* binary_file = fopen(reachable_filename.c_str(), "wb");
        size_t written = 0;

        written = fwrite(&reachable_workfunctions, sizeof(uint64_t), 1, binary_file);
        if (written != 1) {
            PRINT_AND_ABORT("The number of reachable workfunctions was not written correctly.");
        }

        written = fwrite(reachable_wfs_arr, sizeof(workfunction<SIZE>), reachable_workfunctions,
                         binary_file);
        if (written != reachable_workfunctions) {
            PRINT_AND_ABORT("The array of reachable work functions was not written correctly.");
        }

        written = fwrite(adjacent_functions_arr, sizeof(std::array<unsigned int, SIZE>), reachable_workfunctions,
                         binary_file);
        if (written != reachable_workfunctions) {
            PRINT_AND_ABORT("The array of work function adjacencies was not written correctly.");
        }

        written = fwrite(update_costs_arr, sizeof(std::array<short, SIZE>), reachable_workfunctions,
                         binary_file);
        if (written != reachable_workfunctions) {
            PRINT_AND_ABORT("The array of update costs for each workfunction was not written correctly.");
        }

        fclose(binary_file);
    }

    void deserialize_reachable(const std::string& reachable_filename) {
        FILE* binary_file = fopen(reachable_filename.c_str(), "rb");
        size_t read = 0;
        read = fread(&reachable_workfunctions, sizeof(uint64_t), 1, binary_file);
        if (read != 1) {
            PRINT_AND_ABORT("The number of reachable workfunctions was not read correctly.");
        }

        reachable_wfs_arr = new workfunction<SIZE>[reachable_workfunctions];
        update_costs_arr = new std::array<short, SIZE>[reachable_workfunctions];
        adjacent_functions_arr = new std::array<unsigned int, SIZE>[reachable_workfunctions];

        read = fread(reachable_wfs_arr, sizeof(workfunction<SIZE>), reachable_workfunctions,
                         binary_file);
        if (read != reachable_workfunctions) {
            PRINT_AND_ABORT("The array of reachable work functions was not read correctly.");
        }

        read = fread(adjacent_functions_arr, sizeof(std::array<unsigned int, SIZE>), reachable_workfunctions,
                         binary_file);
        if (read != reachable_workfunctions) {
            PRINT_AND_ABORT("The array of work function adjacencies was not read correctly.");
        }

        read = fread(update_costs_arr, sizeof(std::array<short, SIZE>), reachable_workfunctions,
                         binary_file);
        if (read != reachable_workfunctions) {
            PRINT_AND_ABORT("The array of min update costs for each workfunction was not read correctly.");
        }

        fclose(binary_file);
    }


    void fill_hash_to_index() {
        hash_to_index.clear();
        for (int i = 0; i < reachable_workfunctions; i++) {
            hash_to_index[hash(&(reachable_wfs_arr[i]))] = i;
        }
    }

    void initialize_reachable(std::string reachable_filename) {
        if (std::filesystem::exists(reachable_filename)) {
            deserialize_reachable(reachable_filename);
            fill_hash_to_index();
        } else {
            initialize_reachable_from_scratch();
            serialize_reachable(reachable_filename);
        }
    }


    void initialize_reachable_from_scratch() {
        std::unordered_set<uint64_t> reachable_hashes;
        std::vector<algorithm_against_pairwise_opt<SIZE>> reachable_wfs_vec;
        std::vector<std::array<unsigned int, SIZE>> adjacent_functions_vec;
        std::vector<std::array<short, SIZE>> update_costs_vec;

        std::vector<std::array<uint64_t, SIZE>> adjacencies_by_hash;

        algorithm_against_pairwise_opt<SIZE> initial = initial_pwf;
        std::queue<algorithm_against_pairwise_opt<SIZE>> q;
        reachable_hashes.insert(hash(&initial));
        q.push(initial);
        while (!q.empty()) {
            algorithm_against_pairwise_opt<SIZE> front = q.front();
            q.pop();
            hash_to_index[hash(&front)] = reachable_wfs_vec.size();
            reachable_wfs_vec.push_back(front);
            std::array<uint64_t, SIZE> adj;
            std::array<short, SIZE> upd_cost;

            for (short req = 0; req < SIZE; req++) {
                algorithm_against_pairwise_opt<SIZE> new_wf = front;
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
            update_costs_vec.push_back(upd_cost);
        }

        // Create a C-style flat array from the vector. This makes it easier to
        // serialize to a file.

        reachable_workfunctions = reachable_wfs_vec.size();
        reachable_wfs_arr = new workfunction<SIZE>[reachable_workfunctions];
        update_costs_arr = new std::array<short, SIZE>[reachable_workfunctions];
        adjacent_functions_arr = new std::array<unsigned int, SIZE>[reachable_workfunctions];

        for (uint64_t i = 0; i < reachable_workfunctions; i++) {
            reachable_wfs_arr[i] = reachable_wfs_vec[i];
            update_costs_arr[i] = update_costs_vec[i];

            const auto& adj_by_hash = adjacencies_by_hash[i];
            std::array<unsigned int, SIZE> adj_by_index;
            for (int j = 0; j < SIZE; j++) {
                adj_by_index[j] = hash_to_index[adj_by_hash[j]];
            }

            adjacent_functions_arr[i] = adj_by_index;
        }
    }


    void print_reachable(const std::string& filename) {
        FILE* outf = fopen(filename.c_str(), "w");
        for (int i = 0;i < reachable_workfunctions; i++) {
            fprintf(outf, "Reachable pairwise work function of id (index in array) %d:\n", i);
            reachable_wfs_arr[i].print(outf);
        }
        fclose(outf);
    }

};
