#pragma once

#include <algorithm>
#include <cassert>
#include <array>
#include <cinttypes>
#include <cstdint>

#include "common.hpp"
#include "wf/workfunction.hpp"

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

    void swap_inplace(short swap_source) {
        assert(swap_source >= 0 && swap_source <= SIZE-2);
        std::swap(data[swap_source], data[swap_source+1]);
    }

    void move_forward_inplace(short element, short target_pos) {
        short pos = position(element);
        while (pos > target_pos) {
            swap_inplace(pos-1);
            pos--;
        }
    }

    void mtf_inplace(short element) {
        move_forward_inplace(element, 0);
    }

    short position(short element) const {
        for (short i = 0; i < SIZE; i++) {
            if (data[i] == element) {
                return i;
            }
        }
        return -1;
    }


    permutation<SIZE> move_forward_copy(short element, short target_pos) const {
        permutation<SIZE> ret(*this);
        short pos = position(element);
        while (pos > target_pos) {
            ret.swap_inplace(pos-1);
            pos--;
        }
        return ret;
    }

    // Performs move to front.
    permutation<SIZE> mtf_copy(short element) const {
        return move_forward_copy(element, 0);
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

    static permutation<SIZE> perm_from_index_quadratic(uint64_t index) {
        permutation<SIZE> ret;

        std::array<bool, SIZE> placed{};

        for (int i = 0; i < SIZE; i++) {
            short relpos = (short) (index / factorial(SIZE-i-1));
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

            ret.data[i] = candidate;
            placed[candidate] = true;

            index %= factorial(SIZE-i-1);
        }

        return ret;
    }
};
