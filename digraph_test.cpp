#include "digraph.hpp"

int main(void) {
    digraph d1;

    d1.add_vertex();
    d1.add_vertex();
    d1.add_vertex();

    d1.add_edge(0,1, 1.0);
    d1.add_edge(1,2, 1.0);
    d1.add_edge(2,0, 1.0);
    d1.bellman_ford();

    digraph d2;

    d2.add_vertex();
    d2.add_vertex();
    d2.add_vertex();

    d2.add_edge(0,1, 1.0);
    d2.add_edge(1,2, 1.0);
    d2.add_edge(2,0, -2.1);
    d2.bellman_ford();

    return 0;

}