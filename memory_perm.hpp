#pragma once


#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <limits>
#include <cstdio>
#include "common.hpp"
#include "permutations.hpp"

// TODO: use the handy classes from wf/permutation.hpp.


short perm_position(const permutation *p, short element) {
    for (short i = 0; i < LISTSIZE; i++) {
        if ((*p)[i] == element) {
            return i;
        }
    }
    return -1;
}


void perm_swap_inplace(permutation *p, short swap_source) {
    assert(swap_source >= 0 && swap_source <= LISTSIZE-2);
    std::swap((*p)[swap_source], (*p)[swap_source+1]);
}

void perm_move_forward(permutation *p, short element, short target_pos) {
    short pos = perm_position(p, element);
    while (pos > target_pos) {
        perm_swap_inplace(p, pos-1);
        pos--;
    }
}

// Performs move to front.
void perm_mtf(permutation *p, short element)  {
    return perm_move_forward(p, element, 0);
}

class memory_perm {
public:

    static constexpr uint64_t max = factorial(LISTSIZE)-1;

    uint64_t data = 0;

    uint64_t access(int pos) const {
        return perm_from_index_quadratic(data)[pos];
    }

    // Sends the requested element to the front.
    void mtf(int request) {
        auto perm = perm_from_index_quadratic(data);
        if (ALG_DEBUG) {
            fprintf(stderr, "Mem before mtf: ");
            print_permutation(&perm);
        }
        perm_mtf(&perm, request);
        if (ALG_DEBUG) {
            fprintf(stderr, "Mem after mtf: ");
            print_permutation(&perm);
        }

        uint64_t new_index = lexindex_quadratic(&perm);
        data = new_index;
    }

    memory_perm recompute(permutation *alg_relabeling) {
        auto perm = perm_from_index_quadratic(data);
        if(ALG_DEBUG) {
            fprintf(stderr, "Memory before relabeling: ");
            print_permutation(&perm);
            fprintf(stderr, "Relabeling: ");
            print_permutation(alg_relabeling);
        }

        permutation resulting;
        for (int i = 0; i < LISTSIZE; i++) {
            resulting[i] = (*alg_relabeling)[perm[i]];
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
        auto perm = perm_from_index_quadratic(data);
        fprintf(stderr, "[");
        for (int i = 0; i < LISTSIZE; i++) {
            if(i > 0) {
                fprintf(stderr, ", ");
            }
            fprintf(stderr, "%hd" , perm[i]);
        }
        fprintf(stderr, "]");

    }
};
