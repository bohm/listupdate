#pragma once

#include <random>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <filesystem>
#include <iostream>
#include <fstream>

#include "permutation_graph.hpp"
#include "workfunction.hpp"
//#include "parallel-hashmap/parallel_hashmap/phmap.h" // The code requires the parallel-hashmap header-only library.
#include "wf/double_zobrist.hpp"
#include "data_structures/char_flat_set.hpp"
#include "data_structures/file_based_queue.hpp"


template <int SIZE>
class wf_manager {
public:
    permutation_graph<SIZE>& pm;

    std::array<std::array<uint64_t, diameter_bound(SIZE) + 1>, factorial[SIZE]>* zobrist;

    uint64_t reachable_workfunctions = 0;
    workfunction<SIZE>* reachable_wfs_arr = nullptr;
    std::array<unsigned int, SIZE>* adjacent_functions_arr = nullptr;
    std::array<short, SIZE>* min_update_costs_arr = nullptr;

    // std::vector<std::array<unsigned int, SIZE>> adjacent_functions;
    // std::vector<std::array<short, SIZE>> min_update_costs;
    std::unordered_map<uint64_t, unsigned int> hash_to_index;

    wf_manager(permutation_graph<SIZE>& p) : pm(p) {
        zobrist = new std::array<std::array<uint64_t, diameter_bound(SIZE) + 1>, factorial[SIZE]>();
        for (int i = 0; i < factorial[SIZE]; i++) {
            for (int v = 0; v < diameter_bound(SIZE) + 1; v++) {
                (*zobrist)[i][v] = rand_64bit();
            }
        }

        initialize_inversions();
    }

    ~wf_manager() {
        delete zobrist;
    }

    static void initialize_inversions() {
        if (!inversions_ready) {
            invs->vals[0] = 0;
            for (int i = 1; i < factorial[TESTSIZE]; i++) {
                invs->vals[i] = diameter_bound(TESTSIZE);
            }

            dynamic_update(invs);
        }
        inversions_ready = true;
    }

    void flat_update(workfunction<SIZE>* wf, short req) {
        for (int i = 0; i < factorial[SIZE]; i++) {
            wf->vals[i] += pm.all_perms[i].position(req);
        }
    }

    void cut_minimum(workfunction<SIZE>* wf) {
        short m = wf->min();
        for (int i = 0; i < factorial[SIZE]; i++) {
            wf->vals[i] -= m;
        }
    }

    // NEW: our hash now never allows the final byte to be 0000. This way we lose some precision, but we make sure
    // that the hash set only reports actual collisions, not false positives with 0000.

    // This change needs to reflect in the other hash functions too.

    static inline uint64_t avoid_0000(uint64_t hash_candidate) {
        unsigned char last_byte = (unsigned char)hash_candidate;
        if (last_byte != 0) {
            return hash_candidate;
        }
        else {
            return hash_candidate | (uint64_t)1;
        }
    }

    uint64_t hash(workfunction<SIZE>* wf) {
        uint64_t ret = 0;
        for (int i = 0; i < factorial[SIZE]; i++) {
            ret ^= (*zobrist)[i][wf->vals[i]];
        }
        return avoid_0000(ret);
    }

    // We use one permutation as a symmetry of the permutahedron to compute
    // what the hash of the symmetric workfunction to the current one is.
    uint64_t hash_under_right_composition(const workfunction<SIZE>* wf, uint64_t perm_id) const {
        uint64_t h = 0;
        for (int i = 0; i < factorial[SIZE]; i++) {
            h ^= (*zobrist)[pm.quick_compose_right(perm_id, i)][wf->vals[i]];
        }
        return avoid_0000(h);
    }

    // First apply the right composition, the mirror using the inverse permutation from the left.
    // The inverse permutation should be the one with id factorial[SIZE]-1.
    uint64_t hash_under_mirrored_composition(const workfunction<SIZE>* wf, uint64_t perm_id) const {
        uint64_t h = 0;
        for (int i = 0; i < factorial[SIZE]; i++) {
            uint64_t compose_and_mirror = pm.
                quick_compose_left(factorial[SIZE] - 1, pm.quick_compose_right(perm_id, i));
            h ^= (*zobrist)[compose_and_mirror][wf->vals[i]];
        }
        return avoid_0000(h);
    }

    // Two functions used primarily to test the above functions.

    workfunction<SIZE> right_composition(const workfunction<SIZE>* wf, uint64_t perm_id) {
        workfunction<SIZE> ret;
        for (int i = 0; i < factorial[SIZE]; i++) {
            // ret.vals[i] = wf->vals[pm.quick_compose_right(perm_id, i)];
            ret.vals[pm.quick_compose_right(perm_id, i)] = wf->vals[i];
        }
        return ret;
    }

    workfunction<SIZE> mirrored_composition(const workfunction<SIZE>* wf, uint64_t perm_id) {
        workfunction<SIZE> ret;
        for (int i = 0; i < factorial[SIZE]; i++) {
            uint64_t compose_and_mirror = pm.
                quick_compose_left(factorial[SIZE] - 1, pm.quick_compose_right(perm_id, i));
            // ret.vals[i] = wf->vals[compose_and_mirror];
            ret.vals[compose_and_mirror] = wf->vals[i];
        }
        return ret;
    }


    // Static, but requires the pg pointer to be populated.
    static void dynamic_update(workfunction<SIZE>* wf) {
        for (short value = 0; value < diameter_bound(SIZE); value++) {
            for (int i = 0; i < factorial[SIZE]; i++) {
                if (wf->vals[i] == value) {
                    for (uint64_t adj : pg->adjacencies[i]) {
                        wf->vals[adj] = std::min((short)(value + 1), wf->vals[adj]);
                    }
                }
            }
        }
    }

    // Old version.

    uint64_t adjacency_explicit(uint64_t wf_index, short request) {
        workfunction<SIZE> new_wf = reachable_wfs_arr[wf_index];
        flat_update(&new_wf, request);
        cut_minimum(&new_wf);
        dynamic_update(&new_wf);
        // new_wf.validate();
        uint64_t new_wf_index = hash_to_index[hash(&new_wf)];
        return new_wf_index;
    }

    uint64_t adjacency(uint64_t wf_index, short request) {
        return adjacent_functions_arr[wf_index][request];
    }

    short update_cost(uint64_t wf_index, short request) {
        return min_update_costs_arr[wf_index][request];
    }

    void print_all_symmetries(const workfunction<SIZE>* wf) {
        for (unsigned long i = 0; i < factorial[SIZE]; i++) {
            fprintf(stderr, "Composition with permutation %lu:\n", i);
            right_composition(wf, i).print();
            fprintf(stderr, "Composition and mirror with permutation %lu:\n", i);
            mirrored_composition(wf, i).print();
        }
    }

    bool any_symmetry_in_reachable(const workfunction<SIZE>* new_wf,
                                   const char_flat_set& reachable_hashes) {
        for (unsigned long i = 0; i < factorial[SIZE]; i++) {
            uint64_t hash_projection = hash_under_right_composition(new_wf, i);
            workfunction<SIZE> tmp = right_composition(new_wf, i);
            assert(hash_projection == hash(&tmp)); // Debug.
            if (reachable_hashes.contains(hash_projection)) {
                return true;
            }
            uint64_t hash_projection_and_mirror = hash_under_mirrored_composition(new_wf, i);
            tmp = mirrored_composition(new_wf, i);
            assert(hash_projection_and_mirror == hash(&tmp)); // Debug.

            if (reachable_hashes.contains(hash_projection_and_mirror)) {
                return true;
            }
        }

        return false;
    }

    // Counts reachable functions without initializing the full set of reachable functions.
    // Gentler on the memory, useful primarily for estimates.

    uint64_t count_reachable() {
        // std::unordered_set<uint64_t> reachable_hashes;
        char_flat_set reachable_hashes(35);

        // phmap::flat_hash_set<uint64_t> reachable_hashes;
        workfunction<SIZE> initial = *invs;
        std::queue<workfunction<SIZE>> q;
        file_based_queue fbq(std::string("queue.bin"));

        reachable_hashes.insert(hash(&initial));
        fbq.push(&initial);
        // q.push(initial);
        // Debug
        // print_all_symmetries(&initial);

        while (!fbq.empty()) {
            // workfunction<SIZE> front = q.front();
            // q.pop();
            workfunction<SIZE>* front = fbq.pop_nonblocking();
            // front->print(stderr);
            for (short req = 0; req < SIZE; req++) {
                workfunction<SIZE> new_wf = *front;
                flat_update(&new_wf, req);
                cut_minimum(&new_wf);
                dynamic_update(&new_wf);
                uint64_t h = hash(&new_wf);
                // if (!reachable_hashes.contains(h)) {
                if (!any_symmetry_in_reachable(&new_wf, reachable_hashes)) {
                    // fprintf(stderr, "New reachable work function found:\n");
                    // new_wf.print();
                    // fprintf(stderr, "---\n");
                    reachable_hashes.insert(h);
                    fbq.push(&new_wf);
                    if (reachable_hashes.insertions % 1000 == 0) {
                        fprintf(stderr, "Reachable hashes now %lu, queue size %zu.\n", reachable_hashes.insertions,
                                fbq.size());
                    }
                }
            }
            delete front;
        }

        reachable_hashes.report_collisions();
        return reachable_hashes.insertions;
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

        written = fwrite(min_update_costs_arr, sizeof(std::array<short, SIZE>), reachable_workfunctions,
                         binary_file);
        if (written != reachable_workfunctions) {
            PRINT_AND_ABORT("The array of min update costs for each workfunction was not written correctly.");
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
        min_update_costs_arr = new std::array<short, SIZE>[reachable_workfunctions];
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

        read = fread(min_update_costs_arr, sizeof(std::array<short, SIZE>), reachable_workfunctions,
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
        }
        else {
            initialize_reachable_from_scratch();
            serialize_reachable(reachable_filename);
        }
    }

    void initialize_reachable_from_scratch() {
        std::vector<workfunction<SIZE>> reachable_wfs_vec;

        std::unordered_set<uint64_t> reachable_hashes;

        std::vector<std::array<uint64_t, SIZE>> adjacencies_by_hash;
        std::vector<std::array<short, SIZE>> min_update_costs_vec;

        workfunction<SIZE> initial = *invs;
        std::queue<workfunction<SIZE>> q;
        reachable_hashes.insert(hash(&initial));
        q.push(initial);
        while (!q.empty()) {
            workfunction<SIZE> front = q.front();
            q.pop();
            hash_to_index[hash(&front)] = reachable_wfs_vec.size();
            reachable_wfs_vec.push_back(front);
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
            min_update_costs_vec.push_back(upd_cost);
        }

        // Create a C-style flat array from the vector. This makes it easier to
        // serialize to a file.

        reachable_workfunctions = reachable_wfs_vec.size();
        reachable_wfs_arr = new workfunction<SIZE>[reachable_workfunctions];
        min_update_costs_arr = new std::array<short, SIZE>[reachable_workfunctions];
        adjacent_functions_arr = new std::array<unsigned int, SIZE>[reachable_workfunctions];

        for (uint64_t i = 0; i < reachable_workfunctions; i++) {
            reachable_wfs_arr[i] = reachable_wfs_vec[i];
            min_update_costs_arr[i] = min_update_costs_vec[i];

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
        for (int i = 0; i < reachable_workfunctions; i++) {
            fprintf(outf, "Reachable work function of id (index in array) %d:\n", i);
            reachable_wfs_arr[i].print(outf);
        }
        fclose(outf);
    }
};
