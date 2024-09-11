#pragma once


#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <limits>
#include <cstdio>
#include "common.hpp"
#include "old_perm_functions.hpp"
#include "permutation.hpp"
#include "workfunction.hpp"

// A strange combo of a permutation and then a mapping
// that is also numerically a permutation, but we see it as a non-decreasing
// map: 0 -> {0,1,2,3}, 1 -> {1,2,3}, 2 -> {2,3} and so on.

template <int NDMAPSIZE> class nd_map {
public:
    std::array<short, NDMAPSIZE> data;

    nd_map(int index) {
        for (int i = NDMAPSIZE-1; i >= 2; i--) {
            data[NDMAPSIZE-1-i] = index / i;
            assert(index/i <= (NDMAPSIZE-1-i));
            index %= i;
        }
    }
};

template<int DOUBLESIZE> class memory_perm {
public:

    static constexpr uint64_t half_fact = factorial[DOUBLESIZE/2];
    static constexpr uint64_t max = half_fact^2 - 1;

    uint64_t data = 0;

    uint64_t access(int pos) const {
        return permutation<PERMSIZE>::perm_from_index_quadratic(data).data[pos];
    }

    // Sends the requested element to the front.
    void mtf(int request) {
        auto perm = permutation<PERMSIZE>::perm_from_index_quadratic(data);
        perm.mtf_inplace(request);

        if (ALG_DEBUG) {
            fprintf(stderr, "Mem before mtf_copy: ");
            perm.print();
        }
        perm.mtf_inplace(request);

        if (ALG_DEBUG) {
            fprintf(stderr, "Mem after mtf_copy: ");
            perm.print();
        }

        uint64_t new_index = perm.id();
        data = new_index;
    }

    memory_perm recompute(array_as_permutation *alg_relabeling) {
        auto perm = permutation<PERMSIZE>::perm_from_index_quadratic(data);
        if(ALG_DEBUG) {
            fprintf(stderr, "Memory before relabeling: ");
            perm.print();
            fprintf(stderr, "Relabeling: ");
            print_permutation(alg_relabeling);
        }

        permutation<PERMSIZE> resulting{};
        for (int i = 0; i < PERMSIZE; i++) {
            resulting.data[i] = (*alg_relabeling)[perm.data[i]];
        }
        uint64_t new_index = resulting.id();

        if(ALG_DEBUG) {
            fprintf(stderr, "Memory after relabeling: ");
            resulting.print();
        }
        memory_perm m; m.data = new_index;
        return m;
    }

    void full_print() {
        auto perm = permutation<PERMSIZE>::perm_from_index_quadratic(data);
        perm.print();
    }
};
