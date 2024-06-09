#pragma once
#include <cstdint>
#include <cinttypes>
#include <array>
#include <cassert>
#include <algorithm>

constexpr uint64_t factorial(uint64_t n) {
    return n <= 1 ? 1 : (n* factorial(n-1));
}

constexpr short diameter_bound(short n) {
    return (n*(n-1))/2;
}

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
};



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
