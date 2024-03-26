#pragma once

#include <array>
#include "memory.hpp"


constexpr uint64_t factorial (uint64_t n) {
    return n <= 1 ? 1 : (n* factorial(n-1));
}


void recompute_alg_perm(permutation *alg_p, permutation *opt_single_swap) {
    for (int i = 0; i <LISTSIZE; i++) {
        (*alg_p)[i] = (*opt_single_swap)[(*alg_p)[i]];
    }
}

memory recompute_memory(memory m, permutation* alg_relabeling) {
    memory new_mem;
    for (int i = 0; i < LISTSIZE; i++) {
        for (int j = i+1; j < LISTSIZE; j++) {
            uint64_t val = m.access(canonical_order[i*LISTSIZE+j]);
            if (val == 1) {
                // fprintf(stderr, "The pair (%d,%d) needs to be relabeled.\n", i, j);
                int relabeled_i = (*alg_relabeling)[i];
                int relabeled_j = (*alg_relabeling)[j];
                // fprintf(stderr, "The relabeled positions are (%d,%d).\n", relabeled_i, relabeled_j);
                int relabeled_min = std::min(relabeled_i, relabeled_j);
                int relabeled_max = std::max(relabeled_i, relabeled_j);
                new_mem.set_true(canonical_order[relabeled_min*LISTSIZE+relabeled_max]);
            }
        }
    }
    return new_mem;
}

void print_memory_info(memory m) {
    for (int i = 0; i < LISTSIZE; i++) {
        for (int j = i + 1; j < LISTSIZE; j++) {
            uint64_t bit = m.access_sorted_pair(i,j);
            if (bit == 1) {
                fprintf(stderr, "The pair (%d,%d) has been flagged.\n", i, j);
            }
        }
    }
}

constexpr permutation identity() {
    permutation ret {0};
    for (short i = 0; i < LISTSIZE; i++) {
        ret[i] = i;
    }
    return ret;
}


constexpr permutation full_inverse() {
    permutation ret {0};
    for (short i = 0; i < LISTSIZE; i++) {
        ret[i] = LISTSIZE - 1 - i;
    }
    return ret;
}

// Lexicographically next permutation. Returns false if perm was the largest one.
// Knuth's algorithm.
// Implementation from https://www.nayuki.io/page/next-lexicographical-permutation-algorithm.
bool increase(permutation *perm) {
    // Find non-increasing suffix
    size_t i = LISTSIZE - 1;
    while (i > 0 && (*perm)[i - 1] >= (*perm)[i]) {
        i--;
    }

    if (i == 0) {
        return false;
    }

    // Find successor to pivot
    size_t j = LISTSIZE - 1;
    while ((*perm)[j] <= (*perm)[i - 1]) {
        j--;
    }
    short temp = (*perm)[i - 1];
    (*perm)[i - 1] = (*perm)[j];
    (*perm)[j] = temp;

    // Reverse suffix
    j = LISTSIZE - 1;
    while (i < j) {
        temp = (*perm)[i];
        (*perm)[i] = (*perm)[j];
        (*perm)[j] = temp;
        i++;
        j--;
    }
    return true;
}

constexpr permutation IDENTITY = identity();
constexpr permutation FULL_INVERSE = full_inverse();

void iterate_over_permutations(void (*permutation_pointer_function)(permutation *)) {
    permutation iterator = IDENTITY;

    do {
        permutation_pointer_function(&iterator);
    }
    while(increase(&iterator));
}

void iterate_over_memory_and_permutation(void (*perm_and_memory_pointer_function)(permutation *, memory)) {
    permutation iterator = IDENTITY;

    do {
        memory m;
        while (m.data <= max_memory) {
            perm_and_memory_pointer_function(&iterator, m);
            m.data++;
        }
    }
    while(increase(&iterator));

}

void print_permutation(permutation *perm, FILE *f = stderr, bool newline = true) {
    fprintf(f, "(");
    for (int i = 0; i < LISTSIZE; i++) {
        fprintf(f, "%hd", (*perm)[i]);
        if (i<LISTSIZE-1) {
            fprintf(f, ",");
        }
    }
    fprintf(f, ")");

    if(newline) {
        fprintf(f, "\n");
    }
}

void print_permutation_and_memory(permutation *perm, memory m) {
    print_permutation(perm);
    print_memory_info(m);
}

// Quadratic lexicographic order, should be good enough for short arrays.
uint64_t lexindex_quadratic(permutation *perm) {
    uint64_t ret = 0;
    for ( int i = 0; i < LISTSIZE; i ++) {
        uint64_t relative_position = 0;
        if (i >= 1) {
            for (int j = i+1; j < LISTSIZE; j++) {
                if ((*perm)[j] < (*perm)[i] ) {
                    relative_position++;
                }
            }
        } else {
            relative_position = (*perm)[i];
        }

        ret += relative_position*factorial(LISTSIZE-1-i);
    }

    return ret;
}

// Performs a trivial swap.
void swap(permutation *perm, int swap_target) {
    // Allow for no-ops.
    if (swap_target == -1) {
        return;
    }

    std::swap( (*perm)[swap_target], (*perm)[swap_target+1]);
}