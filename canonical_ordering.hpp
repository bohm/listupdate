#pragma once
#include <array>

// using cost_t = long int;
using cost_t = long double;

constexpr unsigned short LISTSIZE = 6;

using permutation = std::array<short, LISTSIZE>;

constexpr std::array<int, LISTSIZE*LISTSIZE> canonical_ordering() {
    std::array<int, LISTSIZE*LISTSIZE> co {0};
    int counter = 0;
    for (int i = 0; i <LISTSIZE; i++) {
        for (int j = i+1; j < LISTSIZE; j++) {
            co[i*LISTSIZE+j] = counter++;
        }
    }
    return co;
}

constexpr std::array<int, LISTSIZE*LISTSIZE> canonical_order = canonical_ordering();

constexpr uint64_t max_memory = (1LLU << (canonical_order[(LISTSIZE-2)*LISTSIZE+(LISTSIZE-1)]+1))-1;