#pragma once

#include <cstdint>
#include <cinttypes>
#include <cassert>
#include <limits>
#include <cstdio>
#include "common.hpp"

constexpr uint64_t ONES = std::numeric_limits<uint64_t>::max();

class memory_pairs {
public:

    static constexpr uint64_t max = max_memory_pairs;

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

    // Higher-order functions for accessing MEMORY with respect to the canonical ordering.

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

    memory_pairs recompute(permutation *alg_relabeling) {
        memory_pairs new_mem;
        for (int i = 0; i < LISTSIZE; i++) {
            for (int j = i+1; j < LISTSIZE; j++) {
                uint64_t val = access(canonical_order[i*LISTSIZE+j]);
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

    void full_print() {
        for (int i = 0; i < LISTSIZE; i++) {
            for (int j = i + 1; j < LISTSIZE; j++) {
                uint64_t bit = access_sorted_pair(i,j);
                if (bit == 1) {
                    fprintf(stderr, "The pair (%d,%d) has been flagged.\n", i, j);
                }
            }
        }
    }
};
