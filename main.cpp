#include "implicit_bellman_ford.hpp"
#include "bellman_ford.hpp"

int main(void) {
    // FILE *f = fopen("cmake-build-debug/graphdump.dot", "w");
    create_graph();

    // long int random_number=1202;
    // auto v = g.get_vert(random_number);
    // assert((long int) v->id == random_number);

    // print_graph(f);
    // fclose(f);

    bellman_ford();
    // fprintf(stderr, "---\n");
    fprintf(stderr, "Implicit computation:\n");
    // bellman_ford_implicit();
    // bellman_ford_implicit_parallel();
    // bellman_ford_compare_exchange();
    // auto [distances_len, distances] = read_distance_array();
    // print_array(distances_len, distances);
    // delete distances;

    return 0;
}
