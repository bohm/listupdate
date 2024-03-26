#include "memory_pairs.hpp"
#include "memory_bitfield.hpp"

int main(void) {

    memory_pairs m;
    m.print_memory();

    m.set_true(0);
    m.print_memory();
    m.set_true(62);
    m.print_memory();
    int j1 = m.access(62);
    fprintf(stderr, "The returned value is %d.\n", j1);
    int j2 = m.access(61);
    fprintf(stderr, "The returned value is %d.\n", j2);
    m.set_false(0);
    m.print_memory();
    m.set_false(1);
    m.print_memory();


    memory_bitfield m2;
    fprintf(stderr, "The possible values of memory are [0, %lu].\n", m2.max);
    return 0;
}