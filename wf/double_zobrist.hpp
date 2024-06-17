#pragma once

#include <cstdint>
#include <unordered_set>
#include <random>
#include "../common.hpp"
#include "workfunction.hpp"

// Mersenne twister.
std::mt19937_64 gen(12345);

uint64_t rand_64bit() {
    uint64_t r = gen();
    return r;
}

class double_hashed_el {
public:
    uint64_t first_part = 0;
    uint64_t second_part = 0;
    bool operator==(const double_hashed_el&) const = default; // since C++20
};


template <int SIZE> class double_zobrist {
    std::array<std::array<uint64_t, diameter_bound(SIZE)+1>, factorial(SIZE)> zobrist_first{};
    std::array<std::array<uint64_t, diameter_bound(SIZE)+1>, factorial(SIZE)> zobrist_second{};

public:

    void init() {
        for (int i = 0; i < factorial(SIZE); i++) {
            for (int v = 0; v < diameter_bound(SIZE)+1; v++) {
                zobrist_first[i][v] = rand_64bit();
                zobrist_second[i][v] = rand_64bit();
            }
        }
    }

    double_hashed_el hash(workfunction<SIZE> *wf) {
            double_hashed_el ret{};
            for (int i = 0; i < factorial(SIZE); i++) {
                ret.first_part ^= zobrist_first[i][wf->vals[i]];
                ret.second_part ^= zobrist_first[i][wf->vals[i]];
            }
            return ret;
    }
};

template<>
struct std::hash<double_hashed_el>
{
    std::size_t operator()(const double_hashed_el& s) const noexcept
    {
        return s.first_part ^ s.second_part;
    }
};