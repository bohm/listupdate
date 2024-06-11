#pragma once

#include <algorithm>
#include <cassert>
#include <array>
#include <cinttypes>
#include <cstdint>

#include "common.hpp"
#include "workfunction.hpp"

template <short SIZE> class permutation {
public:
    std::array<short, SIZE> data;

    uint64_t id() const {
        // Quadratic lexicographic order, should be good enough for short arrays.
        uint64_t ret = 0;
        for ( int i = 0; i < SIZE; i ++) {
            uint64_t relative_position = 0;
                if (i >= 1) {
                    for (int j = i+1; j < SIZE; j++) {
                        if (data[j] < data[i] ) {
                            relative_position++;
                        }
                    }
                } else {
                    relative_position = data[i];
                }

                ret += relative_position*factorial(SIZE-1-i);
            }

            return ret;
        }

    void print(FILE *f = stderr, bool newline = true) {
        fprintf(f, "%" PRIu64 ": (", id());
        for (int i = 0; i < SIZE; i++) {
            fprintf(f, "%hd", data[i]);
            if (i<SIZE - 1) {
                fprintf(f, ",");
            }
        }
        fprintf(f, ")");

        if(newline) {
            fprintf(f, "\n");
        }
    }

    // Performs a trivial swap.
    permutation<SIZE> swap(int swap_source) const {
        // Allow for no-ops.
        // if (swap_source == -1) {
        //     return;
        // }
        assert(swap_source >= 0 && swap_source <= SIZE-2);
        permutation<SIZE> copy = *this;
        std::swap( copy.data[swap_source], copy.data[swap_source+1]);
        return copy;
    }

    short position(short element) const {
        for (short i = 0; i < SIZE; i++) {
            if (data[i] == element) {
                return i;
            }
        }
        return -1;
    }

    short inversions() const {
        return invs.vals[id()];
    }

    short inversions_wrt(const permutation<SIZE> *other) const {
        std::array<short, SIZE> inverse;
        for (int i = 0; i < SIZE; i++) {
            inverse[data[i]] = i;
        }

        permutation<SIZE>copy(*other);
        for (int i = 0; i < SIZE; i++) {
            copy.data[i] = inverse[copy.data[i]];
        }

        return copy.inversions();
    }
};
