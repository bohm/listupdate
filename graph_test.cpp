#include "common.hpp"
#include "permutation_graph.hpp"
#include "graphs/graph_bipartite_mru.hpp"
#include "graphs/bip_mru_potentials.hpp"

int main() {
    pg = new permutation_graph<LISTSIZE>();
    pg->init();
    invs = new workfunction<LISTSIZE>{};
    wf_manager<LISTSIZE>::initialize_inversions();

    gbm.pm = pg;
    gbm.populate_vertices();
    gbm.add_all_outedges();
    // gbm.specific_algorithm_outedges();
    gbm.request_improving_outedges();
    // gbm.print(stderr);

    bip_mru_potentials pots(&gbm);
    pots.compute_equivalence_classes();
    // pots.print_equivalence_class(0);
    // fprintf(stderr, "---\n");
    // pots.print_equivalence_class(1);

    bool anything_updated = true;
    uint64_t iter_count = 0;
    while(anything_updated) {
        //if (iter_count % 10 == 0) {
        fprintf(stderr, "Iteration %" PRIu64 ".\n", iter_count);
        //}
        bool adv_updated = pots.update_adv();
        // bool alg_updated = pots.update_alg();
        bool alg_updated = pots.update_alg_equivalence_classes();
        anything_updated = adv_updated || alg_updated;
        if (gbm.min_adv_potential() >= 1.0) {
            fprintf(stdout, "The min ADV potential is higher than one.\n");
            // gbm.print(stderr);

            return 0;
        }
        iter_count++;
    }

    fprintf(stdout, "The potentials have stabilized with min potential 0.\n");
    // gbm.print(stderr);
    /*
    bool anything_updated = true;
    uint64_t iter_count = 0;


    while(anything_updated) {
        //if (iter_count % 10 == 0) {
        fprintf(stderr, "Iteration %" PRIu64 ".\n", iter_count);
        //}
        bool adv_updated = pots.update_adv();
        bool alg_updated = pots.update_alg();
        anything_updated = adv_updated || alg_updated;
        if (gbm.min_adv_potential() >= 1.0) {
            fprintf(stdout, "The min ADV potential is higher than one.\n");
            return 0;
        }
        iter_count++;
    }

    fprintf(stdout, "The potentials have stabilized with min potential 0.\n");
     */
    return 0;
}