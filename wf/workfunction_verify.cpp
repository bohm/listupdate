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


    std::string reachable_after_decisions_filename = std::string("reachable-after-decisions-") + std::to_string(LISTSIZE)
        + std::string("-ratio-") + std::to_string(RATIO) + std::string(".bin");
    std::string last_three_after_decisions_filename = std::string("last-three-after-") + std::to_string(LISTSIZE)
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
    uint64_t rchbl = wm.reachable_workfunctions;
    fprintf(stderr, "Reachable: %" PRIu64 ".\n", rchbl);

    // The actual deal.
    std::string bin_name{};
    game_graph<TESTSIZE> g(wm, true, bin_name);
    g.build_wfa_minima();
    // Saving last three positions.

    if (std::filesystem::exists(reachable_after_decisions_filename)
        && std::filesystem::exists(last_three_after_decisions_filename)) {
        g.deserialize_reachable_arrays(reachable_after_decisions_filename);
        g.deserialize_decisions(last_three_after_decisions_filename);
        fprintf(stderr, "Solving by using the loaded decision map.\n");
        g.lowerbound_via_last_choices(false);
    } else {
        g.init_last_three();
        if (std::filesystem::exists(last_three_filename)) {
            fprintf(stderr, "Loading last three ADV choices from the binary file.\n");
            g.deserialize_last_three(last_three_filename);
        }

        g.deserialize_reachable_arrays(reachable_vertices_filename);
        fprintf(stderr, "Solving by converting the last three suggestions into decision map.\n");
        g.lowerbound_via_last_choices(true);

    }
    // g.print_top_three_for_reachable();
    // fprintf(stderr, "---\n");
    // g.lowerbound_via_decisions();
    return 0;
}
