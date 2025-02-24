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
    std::string graph_binary_filename = std::string("wfs-graph-") + std::to_string(LISTSIZE)
        + std::string("-ratio-") + std::to_string(RATIO) + std::string(".bin");

    std::string workfunctions_binary_filename = std::string("wfs-reachable-") + std::to_string(LISTSIZE) +
        std::string(".bin");
    std::string last_three_filename = std::string("last-three-maximizers-") + std::to_string(LISTSIZE) +
        std::string(".bin");
    std::string reachable_vertices_filename = std::string("reachable-subgraph-") + std::to_string(LISTSIZE)
        + std::string("-ratio-") + std::to_string(RATIO) + std::string(".bin");

    pg = new permutation_graph<LISTSIZE>();
    pg->init();

    invs = new workfunction<LISTSIZE>{};
    wf_manager<LISTSIZE>::initialize_inversions();
    pg->populate_quick_inversions();
    fprintf(stderr, "Total permutations %zu.\n", pg->all_perms.size());

    wf_manager<TESTSIZE> wm(*pg);

    // invs->print();
    wm.initialize_reachable(workfunctions_binary_filename);
    uint64_t rchbl = wm.reachable_wfs.size();
    fprintf(stderr, "Reachable: %" PRIu64 ".\n", rchbl);

    // The actual deal.
    std::string bin_name{};
    game_graph<TESTSIZE> g(wm, true, bin_name);
    g.build_wfa_minima();
    // Saving last three positions.
    g.init_last_three();
    if (std::filesystem::exists(last_three_filename)) {
        fprintf(stderr, "Loading last three ADV choices from the binary file.\n");
        g.deserialize_last_three(last_three_filename);
    }

    assert(std::filesystem::exists(reachable_vertices_filename));
    g.deserialize_reachable_arrays(reachable_vertices_filename);
    // g.print_top_three_for_reachable();
    // fprintf(stderr, "---\n");
    // g.lowerbound_via_decisions();
    // g.lowerbound_via_last_choices();


    // Debug.
    // g.reset_potentials();
    // fprintf(stderr, "---.\n");
    bool anything_updated = true;
    uint64_t iter_count = 0;
    while(anything_updated) {
        fprintf(stderr, "Iteration %" PRIu64 ".\n", iter_count);
        if (g.reachable_min_adv_potential() <= 0) {
            bool adv_updated = g.reachable_update_adv_only_use_last_three(iter_count);
            bool alg_updated = g.reachable_update_alg_wfa_faster();
            anything_updated = adv_updated || alg_updated;
        }

        if (g.reachable_min_adv_potential() >= 1) {
            fprintf(stdout, "The min ADV potential is higher than one.\n");
            return 0;
        }
        iter_count++;
    }

    fprintf(stdout, "The potentials have stabilized with min potential 0.\n");
    return 0;
}
