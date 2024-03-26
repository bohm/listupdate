//
// Created by bohm on 19.3.24.
//

#include "memory.hpp"

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

    return 0;
}