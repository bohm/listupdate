#include <cstdlib>
#include <cstdio>
#include <unordered_set>
#include <queue>
#include "../permutation_graph.hpp"
#include "../wf_manager.hpp"
#include "../workfunction.hpp"
#include "game_graph.hpp"

int main() {
    std::string workfunctions_filename = std::string("wfs-") + std::to_string(LISTSIZE) + std::string(".log");

    std::string workfunctions_binary_fileanme = std::string("wfs-reachable-") + std::to_string(LISTSIZE) +
        std::string(".bin");
    pg = new permutation_graph<LISTSIZE>();
    pg->init();

    invs = new workfunction<LISTSIZE>{};
    wf_manager<LISTSIZE>::initialize_inversions();
    pg->populate_quick_inversions();
    fprintf(stderr, "Total permutations %zu.\n", pg->all_perms.size());

    wf_manager<TESTSIZE> wm(*pg);

    // invs->print();
    wm.initialize_reachable(workfunctions_binary_fileanme);
    uint64_t rchbl = wm.reachable_wfs.size();
    fprintf(stderr, "Reachable: %" PRIu64 ".\n", rchbl);

    // The actual deal.

    game_graph<TESTSIZE> g(wm, true);
    g.build_wfa_minima();
    bool anything_updated = true;
    uint64_t iter_count = 0;
    while(anything_updated) {
        //if (iter_count % 10 == 0) {
            fprintf(stderr, "Iteration %" PRIu64 ".\n", iter_count);
        //}
        bool adv_updated = g.update_adv();
        // bool alg_updated = g.update_alg();
        bool alg_updated = g.update_alg_stay_or_mtf();
        // bool alg_updated = g.update_alg_single_swap();
        // bool alg_updated = g.update_alg_request_moves_forward();
        // bool alg_updated = g.update_alg_wfa();
        // bool alg_updated = g.update_alg_wfa_faster();
        anything_updated = adv_updated || alg_updated;
        if (g.min_adv_potential() >= 1) {
            fprintf(stdout, "The min ADV potential is higher than one.\n");
            // wm.print_reachable(workfunctions_filename);
            // g.print_potential();

            return 0;
        }
        iter_count++;
    }

    fprintf(stdout, "The potentials have stabilized with min potential 0.\n");
    // wm.print_reachable(workfunctions_filename);
    // g.print_potential();
    return 0;
}
