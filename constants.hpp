#pragma once
#include <string>

#ifndef COMP_RATIO
#error "The float macro constant COMP_RATIO must be passed by the compiler."
#define COMP_RATIO 3.0 // This line is a hack to make G++ spit out only the error above.
#endif

#ifndef TSIZE
#error "The integer list size constant TSIZE must be passed by the compiler."
#define TSIZE 3
#endif


constexpr int TESTSIZE = TSIZE;
constexpr unsigned short LISTSIZE = TSIZE;

constexpr float EPSILON = 0.0001;
constexpr int MULTIPLIER = 100;
constexpr int ALG_MULTIPLIER = 2;
constexpr int ADV_MULTIPLIER = 1;

constexpr long double RATIO = COMP_RATIO;
// constexpr long double RATIO = 2.9;
// constexpr long double RATIO = 3.0;
// constexpr long double RATIO = 3.1;
//constexpr long double RATIO = 3.556;
constexpr float RECENCY_RATIO = 0.6;

// Filenames

std::string workfunctions_log_filename = std::string("wfs-") + std::to_string(LISTSIZE) + std::string(".log");
std::string graph_binary_filename = std::string("wfs-graph-") + std::to_string(LISTSIZE)
    + std::string("-ratio-") + std::to_string(RATIO) + std::string(".bin");
std::string reachable_workfunctions_filename = std::string("wfs-reachable-v2-") + std::to_string(LISTSIZE) +
    std::string(".bin");
std::string last_three_filename = std::string("last-three-maximizers-") + std::to_string(LISTSIZE) +
    std::string(".bin");
std::string reachable_vertices_filename = std::string("reachable-subgraph-") + std::to_string(LISTSIZE)
+ std::string("-ratio-") + std::to_string(RATIO) + std::string(".bin");
std::string reachable_after_decisions_filename = std::string("reachable-after-decisions-") + std::to_string(LISTSIZE)
    + std::string("-ratio-") + std::to_string(RATIO) + std::string(".bin");
std::string last_three_after_decisions_filename = std::string("last-three-after-") + std::to_string(LISTSIZE)
    + std::string("-ratio-") + std::to_string(RATIO) + std::string(".bin");

std::string pairwise_workfunctions_filename = std::string("pwfs-") + std::to_string(LISTSIZE) + std::string(".log");
std::string pairwise_potentials_filename = std::string("pwfs-pots-") + std::to_string(LISTSIZE) + std::string(".log");
std::string pairwise_workfunctions_binary_filename = std::string("pwfs-reachable-v2-") +
    std::to_string(LISTSIZE) + std::string(".bin");