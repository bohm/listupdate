
#include "permutations.hpp"

void print_canonical_order() {
    fprintf(stderr, "LISTSIZE is %d.\n", LISTSIZE);
    for (int i =0; i < LISTSIZE; i++) {
        for (int j = i+1; j < LISTSIZE; j++) {
            fprintf(stderr, "canonical_order(%d,%d) = %d.\n", i,j,canonical_order[i*LISTSIZE+j]);
        }
    }
}
int main(void) {
    print_canonical_order();
    fprintf(stderr, "Memory can be any integer from 0 to %" PRIu64 ".\n", max_memory);
    memory m;
    m.flag_sorted_pair(2,3);
    m.flag_sorted_pair(1,3);
    m.flag_sorted_pair(0,1);
    print_memory_info(m);

    // Suppose ALG goes from 3 1 2 0 after relabeling to 0 2 3 1.
    // The relabel mapping (inverse) is 1 2 3 0. (0 is now called 1, 1 is now called 2, 2 is now called 3, etc.)

    // if the old flagged pairs are (2,3), (1,3) and (0,1), the new flagged pairs should be
    // (3,0),(2,0) and (1,2) -- actually (0,3) and (0,2) and (1,2).

    permutation test_inverse = {1,2,3,0};

    memory m2 = recompute_memory(m,&test_inverse);
    print_memory_info(m2);

    iterate_over_permutations(print_permutation);

    iterate_over_memory_and_permutation(print_permutation_and_memory);
    return 0;
}