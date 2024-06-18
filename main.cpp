// #include "implicit_bellman_ford.hpp"
#include "bellman_ford.hpp"
#include "permutation_graph.hpp"
#include "implicit_bellman_ford.hpp"

int main() {
    // FILE *f = fopen("cmake-build-debug/graphdump.dot", "w");
    pg = new permutation_graph<LISTSIZE>();
    pg->init();
    invs = new workfunction<LISTSIZE>{};
    wf_manager<LISTSIZE>::initialize_inversions();

    /*
    permutation<5> one{{3,1,0,4,2}};
    permutation<5> mru{{0,3,4,2,1}};
    for (int i = 4; i >= 0; i--) {
        permutation<5> moved = one.move_from_position_to_position(4,i);
        int inversions = moved.inversions_wrt(&mru);
        fprintf(stderr, "Inversions %d: ", inversions);
        moved.print();

    }

    memory_perm mem; mem.data = mru.id();
    std::array<short, 5> one_data = one.data;
    alg_single_step_mru_minimize_inv(&one_data, &mem, 2);

    */
    // create_graph();
    // g.dfs_reachability();

    // long int random_number=1202;
    // auto v = g.get_vert(random_number);
    // assert((long int) v->id == random_number);

    // print_graph(f);
    // fclose(f);

    // bellman_ford();
    // fprintf(stderr, "---\n");
    // fprintf(stderr, "Implicit computation:\n");
    // bellman_ford_implicit();
    // bellman_ford_implicit_parallel();
    bellman_ford_compare_exchange();
    // auto [distances_len, distances] = read_distance_array();
    // print_array(distances_len, distances);
    // delete distances;

    return 0;
}
