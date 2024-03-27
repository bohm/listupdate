
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
    fprintf(stderr, "Memory can be any integer from 0 to %" PRIu64 ".\n", max_memory_pairs);
    memory_pairs m;
    m.flag_sorted_pair(2,3);
    m.flag_sorted_pair(1,3);
    m.flag_sorted_pair(0,1);
    m.full_print();


    // Suppose ALG goes from 3 1 2 0 after relabeling to 0 2 3 1.
    // The relabel mapping (inverse) is 1 2 3 0. (0 is now called 1, 1 is now called 2, 2 is now called 3, etc.)

    // if the old flagged pairs are (2,3), (1,3) and (0,1), the new flagged pairs should be
    // (3,0),(2,0) and (1,2) -- actually (0,3) and (0,2) and (1,2).

    permutation test_inverse = {1,2,3,0};

    memory_pairs m2 = m.recompute(&test_inverse);
    m2.full_print();
    fprintf(stderr, "---\n");

    permutation inv_test_p = {3,1,2,0};
    memory_pairs inv_test_m;
    inv_test_m.flag_sorted_pair(2,3);
    inv_test_m.flag_sorted_pair(1,3);
    inv_test_m.flag_sorted_pair(0,1);

    print_permutation_and_memory(&inv_test_p, m);
    permutation opt_p = IDENTITY;
    int opt_position_swap = 1;
    fprintf(stderr, "OPT swaps positions %d and %d.\n", opt_position_swap, opt_position_swap+1);

    swap(&opt_p, opt_position_swap);

    recompute_alg_perm(&inv_test_p, &opt_p);
    memory_pairs inv_test_m2 = inv_test_m.recompute(&opt_p);
    print_permutation_and_memory(&inv_test_p, inv_test_m2);
    fprintf(stderr, "---\n");


    // iterate_over_permutations(print_permutation);
    // iterate_over_memory_and_permutation(print_permutation_and_memory);

    // ---

    std::array<permutation, 4> perms;
    perms[0] = {0, 1, 2, 3};
    perms[1] = {0, 3, 2 ,1};
    perms[2] = {2, 0, 3, 1};
    perms[3] = {3, 2, 1, 0};

    fprintf(stderr, "Lexindex of p1, p2, p3, p4: %" PRIu64 ", %" PRIu64 ", %" PRIu64 ", %" PRIu64 ".\n",
            lexindex_quadratic(&perms[0]), lexindex_quadratic(&perms[1]), lexindex_quadratic(&perms[2]),
            lexindex_quadratic(&perms[3]) );

    permutation swaptest = {2, 0, 3 ,1};

    swap(&swaptest, 1);
    print_permutation(&swaptest);

    fprintf(stderr, "---\n");
    
    for (int i = 0; i < perms.size(); i++ )
    {
        permutation from_index = perm_from_index_quadratic(lexindex_quadratic(&perms[i]));
        print_permutation(&from_index);
    }

    return 0;

}