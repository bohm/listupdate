#include <cstdlib>
#include <cstdio>
#include <unordered_set>
#include <queue>
#include "../perm_manager.hpp"
#include "wf_manager.hpp"
#include "workfunction.hpp"
#include "game_graph.hpp"

int main() {
    perm_manager<TESTSIZE> pm{};
    pm.populate_all_perms();
    fprintf(stderr, "Total permutations %zu.\n", pm.all_perms.size());
    for (auto & all_perm : pm.all_perms) {
        all_perm.print();
    }

    pm.populate_adjacencies();
    pm.print_adjacencies();

    wf_manager<TESTSIZE> wm(pm);

    invs.print();

#if TSIZE == 4
    // Test.
    permutation<TESTSIZE> random{};
    random.data = {0,3,1,2};
    permutation<TESTSIZE> initial{};
    initial.data = {0,1,2,3};
    permutation<TESTSIZE> third{};
    third.data = {3,1,2,0};

    fprintf(stderr, "[0,3,1,2] inversions: %hd.\n", random.inversions());
    fprintf(stderr, "[0,1,2,3] inversions w.r.t. [0,3,1,2]: %hd.\n", initial.inversions_wrt(&random));
    fprintf(stderr, "[0,3,1,2] inversions w.r.t. [0,3,1,2]: %hd.\n", random.inversions_wrt(&random));
    fprintf(stderr, "[3,1,2,0] inversions w.r.t. [0,3,1,2]: %hd.\n", third.inversions_wrt(&random));
#endif

    wm.initialize_reachable();
    uint64_t rchbl = wm.reachable_wfs.size();
    fprintf(stderr, "Reachable: %" PRIu64 ".\n", rchbl);

    // The actual deal.

    game_graph<TESTSIZE> g(wm);
    bool anything_updated = true;
    uint64_t iter_count = 0;
    while(anything_updated) {
        //if (iter_count % 10 == 0) {
            fprintf(stderr, "Iteration %" PRIu64 ".\n", iter_count);
        //}
        bool adv_updated = g.update_adv();
        // bool alg_updated = g.update_alg();
        //bool alg_updated = g.update_alg_mtf();
        // bool alg_updated = g.update_alg_single_swap();
        bool alg_updated = g.update_alg_request_moves_forward();
        anything_updated = adv_updated || alg_updated;
        if (g.min_adv_potential() >= 1) {
            fprintf(stdout, "The min ADV potential is higher than one.\n");
            return 0;
        }
        iter_count++;
    }

    fprintf(stdout, "The potentials have stabilized with min potential 0.\n");
    g.print_potential();
    return 0;
}