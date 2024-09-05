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

class memory_perm {
public:

    static constexpr uint64_t max = factorial[LISTSIZE]-1;

    uint64_t data = 0;

    uint64_t access(int pos) const {
        return permutation<LISTSIZE>::perm_from_index_quadratic(data).data[pos];
    }

    // Sends the requested element to the front.
    void mtf(int request) {
        auto perm = permutation<LISTSIZE>::perm_from_index_quadratic(data);
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
        auto perm = permutation<LISTSIZE>::perm_from_index_quadratic(data);
        if(ALG_DEBUG) {
            fprintf(stderr, "Memory before relabeling: ");
            perm.print();
            fprintf(stderr, "Relabeling: ");
            print_permutation(alg_relabeling);
        }

        permutation<LISTSIZE> resulting{};
        for (int i = 0; i < LISTSIZE; i++) {
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
        auto perm = permutation<LISTSIZE>::perm_from_index_quadratic(data);
        perm.print();
    }
};
