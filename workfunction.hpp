#pragma once

#include <algorithm>
#include <cstdarg>
#include "common.hpp"

// A fprintf(stderr) version of the #error macro.
void PRINT_AND_ABORT(const char *format, ...) {
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stderr, format, argptr);
    va_end(argptr);

    abort();
}


template <int SIZE> class workfunction {
public:
    std::array<short, factorial[SIZE]> vals;

    short min() const {
        return (*std::min_element(vals.begin(), vals.end()));
    }

    short max() const {
        return (*std::max_element(vals.begin(), vals.end()));
    }

    void validate() const {
        for (int i = 0; i < factorial[SIZE]; i++) {
            assert(vals[i] >= 0 && vals[i] <= diameter_bound(SIZE));
        }
    }

    void print(FILE *outf = stderr) const {
        for (int i = 0; i < factorial[SIZE]; i++) {
            fprintf(outf, "wf[%d] = %hd.\n", i, vals[i]);
        }
    }

    static workfunction<SIZE>* serialized_read(FILE* reading_file) {
        int succ_read_shorts = 0;
        auto* ret = new workfunction<SIZE>();
        succ_read_shorts = fread(ret->vals.data(), sizeof(short), factorial[SIZE], reading_file);
        if (succ_read_shorts != factorial[SIZE]) {
            delete ret;
            PRINT_AND_ABORT("Serialized read failed: not enough shorts in the file.\n");
            return nullptr;
        }
        int read_delimeter = 0;
        short delimeter = -1;
        read_delimeter = fread(&delimeter, sizeof(short), 1, reading_file);
        if (delimeter != -1 || read_delimeter != 1) {
            delete ret;
            PRINT_AND_ABORT("Serialized read failed: missing delimeter.\n");
            return nullptr;
        }

        return ret;
    }

    static workfunction<SIZE>* buffer_serialized_read(const short *buf, int starting_pos) {
        auto *ret = new workfunction<SIZE>();
        for (int i = 0; i < factorial[SIZE]; i++) {
            ret->vals[i] = buf[starting_pos+i];
        }
        assert(buf[starting_pos+factorial[SIZE]] == -1);
        return ret;
    }

    void serialized_write(FILE* writing_file) const {
        size_t written = fwrite(vals.data(), sizeof(short), factorial[SIZE], writing_file);
        if (written != factorial[SIZE]) {
            PRINT_AND_ABORT("Serialized write failed: not all shorts got written.\n");
        }
        short del = -1;
        fwrite(&del, sizeof(short), 1, writing_file);
        fflush(writing_file);
    }

    void buffer_serialized_write(short *buf, int starting_pos) const {
        for (int i = 0; i < factorial[SIZE]; i++) {
            buf[starting_pos+i] = vals[i];
        }
        buf[starting_pos + factorial[SIZE]] = -1;
    }
};

workfunction<TESTSIZE> *invs = nullptr;
bool inversions_ready = false;