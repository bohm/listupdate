#include "../common.hpp"
#include "../pairwise_wf_manager.hpp"
#include "../permutation_graph.hpp"
#include "../wf_manager.hpp"
#include "pairwise_game_graph.hpp"

int main(void) {
    std::string workfunctions_filename = std::string("pwfs-") + std::to_string(LISTSIZE) + std::string(".log");
    std::string potentials_filename = std::string("pwfs-pots-") + std::to_string(LISTSIZE) + std::string(".log");

    std::string pairwise_workfunctions_binary_filename = std::string("pwfs-reachable-") +
        std::to_string(LISTSIZE) + std::string(".bin");

    pg = new permutation_graph<LISTSIZE>();
    pg->init();


    invs = new workfunction<LISTSIZE>{};
    wf_manager<LISTSIZE>::initialize_inversions();


    pairwise_wf_manager<LISTSIZE> wfm;

    fprintf(stderr, "Total permutations %zu.\n", pg->all_perms.size());
    for (auto & all_perm : pg->all_perms) {
        all_perm.print();
    }

    wfm.initialize_reachable(pairwise_workfunctions_binary_filename);
    uint64_t rchbl = wfm.reachable_wfs.size();
    fprintf(stderr, "Reachable: %" PRIu64 ".\n", rchbl);


    pairwise_game_graph<TESTSIZE> g(wfm, *pg);
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
        anything_updated = adv_updated || alg_updated;
        if (g.min_adv_potential() >= 1) {
            fprintf(stdout, "The min ADV potential is higher than one.\n");
            // wfm.print_reachable(workfunctions_filename);
            // g.print_potential(potentials_filename);

            return 0;
        }
        iter_count++;
    }

    fprintf(stdout, "The potentials have stabilized with min potential 0.\n");
    // wfm.print_reachable(workfunctions_filename);
    // g.print_potential(potentials_filename);
    return 0;
}