#include "graph.hpp"

int main(void) {
    FILE *f = fopen("cmake-build-debug/graphdump.dot", "w");
    create_graph();

    // long int random_number=1202;
    // auto v = g.get_vert(random_number);
    // assert((long int) v->id == random_number);

    print_graph(f);
    fclose(f);

    bellman_ford();
    return 0;
}