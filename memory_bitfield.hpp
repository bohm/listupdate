#pragma once

#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <limits>
#include <cstdio>
#include "common.hpp"

class memory_bitfield {
public:

    static constexpr uint64_t max = (1LLU << LISTSIZE) -1;

    uint64_t data = 0;

    uint64_t access(int pos) const {
        assert(pos >= 0 && pos < LISTSIZE);
        uint64_t mask = 1LLU << pos;
        uint64_t ret = (data & mask) >> pos;
        assert(ret == 0 || ret == 1);
        return ret;
    }

    void set_true(int pos) {
        assert(pos >= 0 && pos < LISTSIZE);
        uint64_t mask = 1LLU << pos;
        data |= mask;
    }

    void set_false(int pos) {
        assert(pos >= 0 && pos < LISTSIZE);
        uint64_t pre_mask = 1LLU << pos;
        uint64_t mask = ONES ^ pre_mask; // Should flip the pos bit to 0, 1 elsewhere.
        data &= mask;
    }


    memory_bitfield recompute(array_as_permutation *alg_relabeling) {
        memory_bitfield new_mem;
        for (int i = 0; i < LISTSIZE; i++) {
            int relabeled_i = (*alg_relabeling)[i];
            auto val = access(i);
            if (val == 1) {
                new_mem.set_true(relabeled_i);
            }
        }

        return new_mem;
    }

    void full_print() {
        for (int i = 0; i < LISTSIZE; i++) {
                fprintf(stderr, "%" PRIu64 "", access(i));
        }
    }
};
