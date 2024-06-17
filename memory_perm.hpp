#pragma once


#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <limits>
#include <cstdio>
#include "common.hpp"
#include "old_perm_functions.hpp"
#include "permutation.hpp"

// TODO: use the handy classes from wf/permutation.hpp.


short perm_position(const array_as_permutation *p, short element) {
    for (short i = 0; i < LISTSIZE; i++) {
        if ((*p)[i] == element) {
            return i;
        }
    }
    return -1;
}


void perm_swap_inplace(array_as_permutation *p, short swap_source) {
    assert(swap_source >= 0 && swap_source <= LISTSIZE-2);
    std::swap((*p)[swap_source], (*p)[swap_source+1]);
}

void perm_move_forward(array_as_permutation *p, short element, short target_pos) {
    short pos = perm_position(p, element);
    while (pos > target_pos) {
        perm_swap_inplace(p, pos-1);
        pos--;
    }
}

// Performs move to front.
void perm_mtf(array_as_permutation *p, short element)  {
    return perm_move_forward(p, element, 0);
}

class memory_perm {
public:

    static constexpr uint64_t max = factorial(LISTSIZE)-1;

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

        array_as_permutation resulting;
        for (int i = 0; i < LISTSIZE; i++) {
            resulting[i] = (*alg_relabeling)[perm.data[i]];
        }
        uint64_t new_index = lexindex_quadratic(&resulting);
        if(ALG_DEBUG) {
            fprintf(stderr, "Memory after relabeling: ");
            print_permutation(&resulting);
        }
        memory_perm m; m.data = new_index;
        return m;
    }

    void full_print() {
        auto perm = permutation<LISTSIZE>::perm_from_index_quadratic(data);
        perm.print();
    }
};
