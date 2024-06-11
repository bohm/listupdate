#pragma once
#include <cstdint>

constexpr uint64_t factorial(uint64_t n) {
    return n <= 1 ? 1 : (n* factorial(n-1));
}

constexpr short diameter_bound(short n) {
    return (n*(n-1))/2;
}

#define TSIZE 3

constexpr bool GRAPH_DEBUG = false;
constexpr int TESTSIZE = TSIZE;
constexpr short RATIO = 3;