//
// Created by boehm on 6/14/24.
//

#include <cstdio>
#include "../common.hpp"
#include "../permutation_graph.hpp"
#include "../wf_manager.hpp"

int main()
{

    pg = new permutation_graph<LISTSIZE>();
    pg->init();
    invs = new workfunction<LISTSIZE>{};
    wf_manager<LISTSIZE>::initialize_inversions();
    permutation_graph<TESTSIZE> pm{};
    pm.populate_all_perms();
    fprintf(stderr, "Total permutations %zu.\n", pm.all_perms.size());
    // for (unsigned long i = 0; i < factorial[LISTSIZE]; i++) {
    //    pm.all_perms[i].print();
    // }
    pm.populate_adjacencies();
    pm.populate_composition();

    wf_manager<LISTSIZE> wm(pm);

    /*workfunction<LISTSIZE> pseudo_wf;
    for (unsigned long u = 0; u < factorial[LISTSIZE]; u++) {
        pseudo_wf.vals[u] = u;
    }
    wm.print_all_symmetries(&pseudo_wf);
     */
    uint64_t rchbl = wm.count_reachable();
    fprintf(stderr, "Reachable work functions: %" PRIu64 ".\n", rchbl);
    return 0;
}
