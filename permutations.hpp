#pragma once

#include <array>
#include "common.hpp"

void recompute_alg_perm(permutation *alg_p, permutation *opt_single_swap) {
    for (int i = 0; i <LISTSIZE; i++) {
        (*alg_p)[i] = (*opt_single_swap)[(*alg_p)[i]];
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

permutation perm_from_index_quadratic(uint64_t index) {
    permutation ret;

    std::array<bool, LISTSIZE> placed{};

    for (int i = 0; i < LISTSIZE; i++) {
        short relpos = (short) (index / factorial(LISTSIZE-i-1));
        // fprintf(stderr, "Iteration %d: Computer relpos %hd.\n", i, relpos);
        short candidate = 0;
        // Compute the i-th unplaced number.
        while (true) {
            if (relpos == 0 && !placed[candidate]) {
                break;
            }

            if (!placed[candidate]) {
                relpos--;
            }
            candidate++;
        }

        // assert(candidate >= 0 && candidate < LISTSIZE && !placed[candidate]);
        // fprintf(stderr, "Iteration %d: Placing digit %hd.\n", i, candidate);

        ret[i] = candidate;
        placed[candidate] = true;

        index %= factorial(LISTSIZE-i-1);
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

permutation inverse(const permutation& p) {
    permutation ret;
    for (short i = 0; i < LISTSIZE; i++) {
        ret[p[i]] = i;
    }
    return ret;
}