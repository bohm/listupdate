#pragma once

#include "memory_pairs.hpp"
#include "memory_bitfield.hpp"
#include "memory_perm.hpp"
#include "old_perm_functions.hpp"
#include "iteration_over_memory.hpp"
#include "permutation.hpp"
#include "wf_manager.hpp"


int alg_single_step_bitfield(array_as_permutation *perm, memory_bitfield *mem, unsigned short presented_item) {
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

    if (FRONT_ACCESS_COSTS_ONE) {
        alg_cost += 1;
    }

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
int alg_single_step_stars(array_as_permutation *perm, memory_pairs *mem, unsigned short presented_item) {
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



// Returns ALG's cost and may edit both permutation and memory.
int alg_single_step_xoror(array_as_permutation *perm, memory_pairs *mem, unsigned short presented_item) {
    int alg_cost = 0;
    int item_pos = 0;
    int target_item_pos = -1;
    for (; item_pos < LISTSIZE; item_pos++) {
        if ((*perm)[item_pos] == presented_item) {
            // if there was no item with bit 0 before presented_item, we will not move.
            if (target_item_pos == -1) {
                target_item_pos = item_pos;
            }
            break;
        } else {
            // Set target_item_pos to the first position from the left which has a bit value 0.
            if (target_item_pos == -1) {
                uint64_t bit = mem->access_pair((*perm)[item_pos], presented_item);
                if (bit == 0) {
                    target_item_pos = item_pos;
                }
            }
        }
    }

    alg_cost += item_pos;

    // bits before j get XORed, bits after j get set to 1 (ORed).
    for (int j = 0; j < LISTSIZE; j++) {
        if (j < item_pos) {
            uint64_t bit = mem->access_pair((*perm)[item_pos], (*perm)[j]);
            if (bit == 1) {
                mem->clear_unsorted_pair((*perm)[item_pos], (*perm)[j]);
            } else {
                mem->flag_unsorted_pair((*perm)[item_pos], (*perm)[j]);
            }
        }

        if (j == item_pos) {
            continue;
        }

        if (j > item_pos) {
            mem->flag_unsorted_pair((*perm)[item_pos], (*perm)[j]);
        }
    }

    // Swap presented_item from item_pos to the left to target_item_pos.

    while(item_pos - target_item_pos > 0) {
        swap(perm, item_pos-1);
        alg_cost++;
        item_pos--;
    }

    return alg_cost;
}



int alg_single_step_mru(array_as_permutation *perm, memory_perm *mem, unsigned short presented_item) {
    if (ALG_DEBUG) {
        fprintf(stderr, "ALG seeking item %d with state: ", presented_item);
        fprintf(stderr, "Memory index %lu.\n", mem->data);
        print_permutation_and_memory<memory_perm>(perm, *mem);
    }

    int alg_cost = 0;
    int item_pos = 0;
    uint64_t flag_cnt = 0;
    for (; item_pos < LISTSIZE; item_pos++) {
        if ((*perm)[item_pos] == presented_item) {
            break;
        }
    }

    alg_cost += item_pos;
    array_as_permutation explicit_memory = perm_from_index_quadratic(mem->data);
    array_as_permutation explicit_memory_inverse = inverse(explicit_memory);

    while (item_pos >= 1 && explicit_memory_inverse[(*perm)[item_pos-1]] > explicit_memory_inverse[presented_item]) {
        swap(perm, item_pos-1);
        alg_cost++;
        item_pos--;
    }

    mem->mtf(presented_item);

    if (ALG_DEBUG) {
        fprintf(stderr, "ALG's cost: %d.\n", alg_cost);
    }
    return alg_cost;
}


int alg_single_step_mru_minimize_inv(array_as_permutation *perm, memory_perm *mem, unsigned short presented_item) {
    if (ALG_DEBUG) {
        fprintf(stderr, "ALG seeking item %d with state: ", presented_item);
        fprintf(stderr, "Memory index %lu.\n", mem->data);
        print_permutation_and_memory<memory_perm>(perm, *mem);
    }

    int alg_cost = 0;
    int item_pos = 0;
    uint64_t flag_cnt = 0;
    for (; item_pos < LISTSIZE; item_pos++) {
        if ((*perm)[item_pos] == presented_item) {
            break;
        }
    }

    alg_cost += item_pos;
    // Compute the position of presented_item, include the search for it in alg_cost.

    permutation<LISTSIZE> mru_memory = permutation<LISTSIZE>::perm_from_index_quadratic(mem->data);

    // Compute the minimum inversion count.
    int min_inversions_val = std::numeric_limits<int>::max();
    int min_inversions_target = item_pos;

    permutation<LISTSIZE> perm_object(*perm);
    for (int move_target = item_pos; move_target >= 0; move_target-- ) {
        auto perm_after_move = perm_object.move_from_position_to_position(item_pos, move_target);
        int inversions = perm_after_move.inversions_wrt(&mru_memory);

        if (ALG_DEBUG) {
            fprintf(stderr, "Candidate for a move (%d inversions): ", inversions);
            perm_after_move.print();
        }
        if (inversions <= min_inversions_val) {
            if(ALG_DEBUG) {
                fprintf(stderr, "The number of inversions %d is lower than previous min of %d.\n",
                        inversions, min_inversions_val);
            }
            min_inversions_val = inversions;
            min_inversions_target = move_target;
        }
    }

    while (item_pos >= 1 && item_pos > min_inversions_target) {
        swap(perm, item_pos-1);
        alg_cost++;
        item_pos--;
    }

    mem->mtf(presented_item);

    if (ALG_DEBUG) {
        fprintf(stderr, "ALG's cost: %d.\n", alg_cost);
    }
    return alg_cost;
}
