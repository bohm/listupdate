#pragma once

#include <algorithm>
#include "common.hpp"

template <int SIZE> class workfunction {
public:
    std::array<short, factorial[SIZE]> vals;

    short min() const {
        return (*std::min_element(vals.begin(), vals.end()));
    }

    short max() const {
        return (*std::max_element(vals.begin(), vals.end()));
    }

    void validate() const {
        for (int i = 0; i < factorial[SIZE]; i++) {
            assert(vals[i] >= 0 && vals[i] <= diameter_bound(SIZE));
        }
    }

    void print() const {
        for (int i = 0; i < factorial[SIZE]; i++) {
            fprintf(stderr, "wf[%d] = %hd.\n", i, vals[i]);
        }
    }
};

workfunction<TESTSIZE> *invs = nullptr;