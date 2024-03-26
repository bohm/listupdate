#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <limits>
#include "memory.hpp"
#include "permutations.hpp"





std::pair<permutation*, uint64_t> edge_to_optstate(algorithm_vertex *v, permutation* new_opt_perm) {

}

void compute_adv_outedges(adversary_vertex *adv_v) {

    for (short item = 0; item < LISTSIZE; item++) {

    }
}
void compute_alg_outedges(algorithm_vertex *alg_v) {
    permutation opt_new_perm = IDENTITY;
    do {
        permutation inverse;
        for (short i = 0; i < LISTSIZE; i++) {
            inverse[opt_new_perm[i]] = i;
        }
        permutation alg_new_perm;
        for (short i = 0; i < LISTSIZE; i++) {
            alg_new_perm[i] = inverse[alg_v->ar[i]];
        }

        memory alg_new_mem = recompute_memory(alg_v->memory, &inverse);

        // Find adv_vertex (memory, permutation);

        // Compute edge cost.

        // Create the edge and add it to the list of edges with the proper cost.

    } while(increase(&opt_new_perm));
}

short inversion_count_quad(permutation* perm) {
    short inversions = 0;
    for (int i = 0; i < LISTSIZE; i++) {
        for (int j = i+1; j < LISTSIZE; j++) {
            if (perm[j] < perm[i]) {
                inversions++;
            }
        }
    }

    return inversions;
}

int algorithmic_step(adversary_vertex *origin, short item) {
    int alg_cost = origin->position(item);


}
template <int LISTSIZE> void compute_one_edge(adversary_vertex<LISTSIZE> *origin, short item,
        std::array<short, LISTSIZE> target_opt_perm) {
    int cost = 0, alg_cost = 0;


}
int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
