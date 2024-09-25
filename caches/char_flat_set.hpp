#pragma once

#include <cinttypes>
#include <cstring>
#include <cassert>
#include <string>

// with some memory allocated dynamically, we also need a dynamic logpart
uint64_t quicklog(uint64_t x) {
    uint64_t ret = 0;
    while (x >>= 1) {
        ret++;
    }
    return ret;
}

inline uint64_t logpart(uint64_t x, int log) {
    return x >> (64 - log);
}

// Powering 2^X.
uint64_t two_to(uint64_t x) {
    return 1LLU << x;
}

// A math routine computing the largest power of two less than x.
uint64_t power_of_two_below(uint64_t x) {
    return two_to(quicklog(x));
}


class char_flat_set
{
public:
    unsigned char *ht;
    uint64_t htsize;
    int logsize;

    uint64_t collisions = 0;
    uint64_t insertions = 0;

     char_flat_set(uint64_t logbytes, std::string descriptor = "") {
        assert(logbytes <= 64);

        uint64_t bytes = two_to(logbytes);
        const uint64_t megabyte = 1024 * 1024;

        htsize = power_of_two_below(bytes / sizeof(unsigned char));
        logsize = quicklog(htsize);
        fprintf(stderr, "Given %lu logbytes (%lu MBs) and el. size %zu, creating %s state cache (64-bit hashes) to %lu els (logsize %d).\n",
                logbytes, bytes / megabyte, sizeof(unsigned char), descriptor.c_str(), htsize, logsize);


        ht = new unsigned char[htsize];
        memset(ht, 0, htsize);
    }

    ~char_flat_set() {
        delete ht;
    }

    uint64_t trim(uint64_t ha) const {
        return logpart(ha, logsize);
    }

    // As long as types are unsigned, this should be the last byte of the uint64_t.
    unsigned char lastchar(uint64_t hash) const {
        return (unsigned char) hash;
    };

    bool contains(uint64_t hash) const {
        unsigned char l = lastchar(hash);
        uint64_t pos = trim(hash);
        return l == ht[pos];

    }

    void insert(uint64_t hash) {
        insertions++;
        unsigned char l = lastchar(hash);
        uint64_t pos = trim(hash);
        if (ht[pos] != 0) {
            collisions++;
        }
        ht[pos] = l;
    }

    uint64_t size() const {
        return htsize;
    }

    void report_collisions() const {
        fprintf(stderr, "Total collisions: %" PRIu64 ".\n", collisions);
    }
};