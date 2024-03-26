#pragma once

#include "memory_pairs.hpp"
#include "permutations.hpp"
#include "memory_bitfield.hpp"


constexpr bool ALG_DEBUG = false;


int alg_single_step_bitfield(permutation *perm, memory_bitfield *mem, unsigned short presented_item) {
    if (ALG_DEBUG) {
        fprintf(stderr, "ALG seeking item %d with state: ", presented_item);
        print_permutation_and_memory<memory_bitfield>(perm, *mem);
    }

    int alg_cost = 0;
    int item_pos = 0;
    uint64_t flag_cnt = 0;
    for (; item_pos < LISTSIZE; item_pos++) {
        if ((*perm)[item_pos] == presented_item) {
            break;
        }
    }

    uint64_t position_bit = mem->access(presented_item);
    if (ALG_DEBUG) {
        fprintf(stderr, "ALG found the item at position %d, and bit value is %lu.\n", item_pos,
                position_bit);
    }

    alg_cost += item_pos;
    if (position_bit == 0) {
        if (item_pos != 0) {
            mem->set_true(presented_item);
        }
        if (ALG_DEBUG && item_pos != 0) {
            fprintf(stderr, "ALG will not swap, the bit was not set.\n");
        } else if (ALG_DEBUG) {
            fprintf(stderr, "ALG has %d in front, it will not swap.\n", presented_item);
        }
    } else {
        if (ALG_DEBUG) {
            fprintf(stderr, "ALG will swap.\n");
        }

        mem->set_false(presented_item);

        while(item_pos > 0) {
            swap(perm, item_pos-1);
            alg_cost++;
            item_pos--;
        }
    }

    if (ALG_DEBUG) {
        fprintf(stderr, "ALG's cost: %d.\n", alg_cost);
    }
    return alg_cost;
}


// Returns ALG's cost and may edit both permutation and memory.
int alg_single_step_original(permutation *perm, memory_pairs *mem, unsigned short presented_item) {
    if (ALG_DEBUG) {
        fprintf(stderr, "ALG seeking item %d with state: ", presented_item);
        print_permutation_and_memory(perm, *mem);
    }
    int alg_cost = 0;
    int item_pos = 0;
    uint64_t flag_cnt = 0;
    for (; item_pos < LISTSIZE; item_pos++) {
        if ((*perm)[item_pos] == presented_item) {
            break;
        } else {
            flag_cnt += mem->access_pair(presented_item, (*perm)[item_pos]);
        }
    }

    if (ALG_DEBUG) {
        fprintf(stderr, "ALG found the item at position %d, with %" PRIu64 " flags on the way.\n", item_pos, flag_cnt);
    }

    uint64_t threshold = (item_pos+1)/2;
    alg_cost += item_pos;
    if (item_pos == 0 || flag_cnt < threshold) {
        // Before returning the alg cost, flag all pairs on the way as true.
        for (int i = 0; i < item_pos; i++) {
            mem->flag_unsorted_pair(presented_item, (*perm)[i]);
        }

        if (ALG_DEBUG && item_pos != 0) {
            fprintf(stderr, "ALG has seen %lu/%lu flags, it will set new flags and not swap.\n", flag_cnt, threshold);
        } else if (ALG_DEBUG) {
            fprintf(stderr, "ALG has %d in front, will not set new flags or swap.\n", presented_item);
        }

    } else {
        if (ALG_DEBUG) {
            fprintf(stderr, "ALG has seen %lu/%lu flags, it swap to the front and clear flags.\n", flag_cnt, threshold);
        }

        // Flag all pairs on the way as false, but then move the item to the front.
        for (int i = 0; i < item_pos; i++) {
            mem->clear_unsorted_pair(presented_item, (*perm)[i]);
        }

        while(item_pos > 0) {
            swap(perm, item_pos-1);
            alg_cost++;
            item_pos--;
        }

    }

    if (ALG_DEBUG) {
        fprintf(stderr, "ALG's cost: %d.\n", alg_cost);
    }
    return alg_cost;
}