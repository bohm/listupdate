#pragma once

#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <limits>
#include <cstdio>
#include "canonical_ordering.hpp"

constexpr uint64_t ONES = std::numeric_limits<uint64_t>::max();

class memory {
public:
    uint64_t data = 0;

    uint64_t access(int pos) const {
        assert(pos >= 0 && pos < 64);
        uint64_t mask = 1LLU << pos;
        uint64_t ret = (data & mask) >> pos;
        assert(ret == 0 || ret == 1);
        return ret;
    }

    void set_true(int pos) {
        assert(pos >= 0 && pos < 64);
        uint64_t mask = 1LLU << pos;
        data |= mask;
    }

    void set_false(int pos) {
        assert(pos >= 0 && pos < 64);
        uint64_t pre_mask = 1LLU << pos;
        uint64_t mask = ONES ^ pre_mask; // Should flip the pos bit to 0, 1 esewhere.
        data &= mask;
    }

    void print_memory(bool newline = true) {
        for (int i = 0; i < 64; i++) {
            fprintf(stderr, "%" PRIu64 "", access(i));
        }

        if(newline) {
            fprintf(stderr, "\n");
        }
    }

    // Higher-order functions for accessing memory with respect to the canonical ordering.

    uint64_t access_sorted_pair(int i, int j) const {
        return access(canonical_order[i*LISTSIZE+j]);
    }

    uint64_t access_pair(int i, int j) const {
        return access_sorted_pair(std::min(i,j), std::max(i,j));
    }

    void flag_sorted_pair(int i, int j) {
        set_true(canonical_order[i*LISTSIZE+j]);
    }

    void flag_unsorted_pair(int i, int j) {
        flag_sorted_pair(std::min(i,j), std::max(i,j));
    }

    void clear_sorted_pair(int i, int j) {
        set_false(canonical_order[i*LISTSIZE+j]);
    }

    void clear_unsorted_pair(int i, int j) {
        clear_sorted_pair(std::min(i,j), std::max(i,j));
    }
};

