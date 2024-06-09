#pragma once

#include <random>
#include "perm_manager.hpp"

template <int SIZE> class workfunction {
public:
    std::array<short, factorial(SIZE)> vals;

    short min() {
        return (*std::min_element(vals.begin(), vals.end()));
    }

    void print() const {
        for (int i = 0; i < factorial(SIZE); i++) {
            fprintf(stderr, "wf[%d] = %hd.\n", i, vals[i]);
        }
    }
};


// Mersenne twister.
std::mt19937_64 gen(12345);

uint64_t rand_64bit() {
    uint64_t r = gen();
    return r;
}

template <int SIZE> class wf_manager {
public:
    perm_manager<SIZE>& pm;

    std::array<std::array<uint64_t, diameter_bound(SIZE)+1>, factorial(SIZE)> zobrist;

    wf_manager(perm_manager<SIZE> &p) : pm(p) {

        for (int i = 0; i < factorial(SIZE); i++) {
            for (int v = 0; v < diameter_bound(SIZE)+1; v++) {
                zobrist[i][v] = rand_64bit();
            }
        }
    }

    void flat_update(workfunction<SIZE> *wf, short req) {
        for (int i = 0; i < factorial(SIZE); i++) {
            wf->vals[i] += pm.all_perms[i].position(req);
        }
    }

    void cut_minimum(workfunction<SIZE> *wf) {
        short m = wf->min();
        for (int i = 0; i < factorial(SIZE); i++) {
            wf->vals[i] -= m;
        }
    }

    uint64_t hash(workfunction<SIZE> *wf) {
        uint64_t ret = 0;
        for (int i = 0; i < factorial(SIZE); i++) {
            ret ^= zobrist[i][wf->vals[i]];
        }
        return ret;
    }

    void dynamic_update(workfunction<SIZE> *wf) {
        for (short value = 0; value < diameter_bound(SIZE); value++) {
            for (int i = 0; i < factorial(SIZE); i++) {
                if (wf->vals[i] == value) {
                    for (uint64_t adj: pm.adjacencies[i]) {
                        wf->vals[adj] = std::min((short) (value+1), wf->vals[adj]);
                    }
                }
            }
        }
    }
};
