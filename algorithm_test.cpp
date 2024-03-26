#include "algorithm.hpp"

int main(void) {

    permutation start = {2, 0, 3, 1};
    memory m;

    m.flag_sorted_pair(2,3);
    m.flag_sorted_pair(1,3);
    m.flag_sorted_pair(0,1);
    print_memory_info(m);

    fprintf(stderr, "---\n");
    alg_single_step(&start, &m, 2);
    fprintf(stderr, "---\n");
    alg_single_step(&start, &m, 3);
    fprintf(stderr, "---\n");

    alg_single_step(&start, &m, 1);
    fprintf(stderr, "---\n");
    print_permutation_and_memory(&start, m);
}