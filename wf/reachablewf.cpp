//
// Created by boehm on 6/14/24.
//

#include <cstdio>
#include "../common.hpp"
#include "../permutation_graph.hpp"
#include "../wf_manager.hpp"

int main()
{

    permutation_graph<TESTSIZE> pm{};
    pm.populate_all_perms();
    fprintf(stderr, "Total permutations %zu.\n", pm.all_perms.size());
    pm.populate_adjacencies();

    wf_manager<TESTSIZE> wm(pm);
    uint64_t rchbl = wm.count_reachable();
    fprintf(stderr, "Reachable: %" PRIu64 ".\n", rchbl);
    return 0;
}