#pragma once
#include "old_perm_functions.hpp"
#include "memory_bitfield.hpp"
#include "memory_pairs.hpp"
#include "memory_perm.hpp"

void iterate_over_memory_and_permutation(void (*perm_and_memory_pointer_function)(array_as_permutation *, MEMORY)) {
    array_as_permutation iterator = IDENTITY;

    do {
        MEMORY m;
        while (m.data <= MEMORY::max) {
            perm_and_memory_pointer_function(&iterator, m);
            m.data++;
        }
    }
    while(increase(&iterator));

}


template <class mem> void print_permutation_and_memory(array_as_permutation *perm, mem m) {
    print_permutation(perm);
    m.full_print();
}
