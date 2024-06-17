#pragma once
#include <cstdint>
#include <cinttypes>
#include <array>
#include <cassert>
#include <algorithm>

#include "common.hpp"
#include "permutation.hpp"


template <short SIZE> class perm_manager
{
public:
    constexpr permutation<SIZE> identity() {
        permutation<SIZE> ret;
        for (short i = 0; i < SIZE; i++) {
            ret.data[i] = i;
        }
        return ret;
    }


    constexpr permutation<SIZE> full_inverse() {
        permutation<SIZE> ret;
        for (short i = 0; i < SIZE; i++) {
            ret[i] = SIZE - 1 - i;
        }
        return ret;
    }

    std::array<permutation<SIZE>, factorial(SIZE)> all_perms;
    std::array<std::array<uint64_t, SIZE-1>, factorial(SIZE)> adjacencies;


    // Lexicographically next permutation. Returns false if perm was the largest one.
    // Knuth's algorithm.
    // Implementation from https://www.nayuki.io/page/next-lexicographical-permutation-algorithm.
    static bool increase(permutation<SIZE> *perm) {
        // Find non-increasing suffix
        size_t i = SIZE - 1;
        while (i > 0 && perm->data[i - 1] >= perm->data[i]) {
            i--;
        }

        if (i == 0) {
            return false;
        }

        // Find successor to pivot
        size_t j = SIZE - 1;
        while (perm->data[j] <= perm->data[i - 1]) {
            j--;
        }
        short temp = perm->data[i - 1];
        perm->data[i - 1] = perm->data[j];
        perm->data[j] = temp;

        // Reverse suffix
        j = SIZE - 1;
        while (i < j) {
            temp = perm->data[i];
            perm->data[i] = perm->data[j];
            perm->data[j] = temp;
            i++;
            j--;
        }
        return true;
    }


    void populate_all_perms() {
        permutation<SIZE> iterator = identity();
        unsigned int i = 0;
        do {
            all_perms[i] = iterator;
            i++;
        }
        while(increase(&iterator));
    }

    void populate_adjacencies() {
        for (int i = 0; i < all_perms.size(); i++) {
            for (int swap = 0; swap < SIZE-1; swap++) {
                uint64_t swapped_id =  all_perms[i].swap(swap).id();
                adjacencies[i][swap] = swapped_id;
            }
        }
    }

    void print_adjacencies() const {
        for (int i = 0; i < adjacencies.size(); i++) {
            fprintf(stderr, "%d: [", i);
            for (int j = 0; j < SIZE-1; j++) {
                if (j > 0) {
                    fprintf(stderr, ",");
                }
                fprintf(stderr, "%" PRIu64, adjacencies[i][j]);
            }
            fprintf(stderr, "]\n");
        }
    }
};
