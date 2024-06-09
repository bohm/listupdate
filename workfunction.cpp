#include <cstdlib>
#include <cstdio>
#include <unordered_set>
#include <queue>
#include "perm_manager.hpp"
#include "wf_manager.hpp"

constexpr int TESTSIZE = 3;

template <short SIZE> uint64_t reachable_work_functions(wf_manager<SIZE> &wm) {
    uint64_t count = 0;
    workfunction<SIZE> initial{};
    initial.vals[0] = 0;
    for (int i = 1; i < factorial(TESTSIZE); i++) {
        initial.vals[i] = diameter_bound(TESTSIZE);
    }

    wm.dynamic_update(&initial);

    std::unordered_set<uint64_t> visited;
    std::queue<workfunction<SIZE>> q;
    q.push(initial);
    while (!q.empty()) {
        workfunction<SIZE> front = q.front();
        q.pop();
        visited.insert(wm.hash(&front));
        count++;
        for (short req = 0; req < SIZE; req++) {
            workfunction<SIZE> new_wf = front;
            wm.flat_update(&new_wf, req);
            wm.cut_minimum(&new_wf);
            wm.dynamic_update(&new_wf);
            uint64_t hash = wm.hash(&new_wf);
            if (!visited.contains(hash)) {
                q.push(new_wf);
            }
        }
    }
    return count;
}

int main() {
    perm_manager<TESTSIZE> pm{};
    pm.populate_all_perms();
    fprintf(stderr, "Total permutations %zu.\n", pm.all_perms.size());
    for (auto & all_perm : pm.all_perms) {
        all_perm.print();
    }

    pm.populate_adjacencies();
    pm.print_adjacencies();

    wf_manager<TESTSIZE> wm(pm);

    workfunction<TESTSIZE> initial{};
    initial.vals[0] = 0;
    for (int i = 1; i < factorial(TESTSIZE); i++) {
        initial.vals[i] = diameter_bound(TESTSIZE);
    }

    wm.dynamic_update(&initial);
    initial.print();

    uint64_t rchbl = reachable_work_functions<TESTSIZE>(wm);
    fprintf(stderr, "Reachable: %" PRIu64 ".\n", rchbl);
    return 0;
}