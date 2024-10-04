#pragma once

#include <cstdio>
#include <string>

#include "../common.hpp"
#include "../workfunction.hpp"

class file_based_queue {
public:
    static constexpr int BUF_WORKFUNCTIONS = 10000;
    static constexpr unsigned int BUF_CAPACITY = (factorial[LISTSIZE]+1)*BUF_WORKFUNCTIONS;
    FILE *reading = nullptr;
    FILE *appending = nullptr;
    int seeking_position = 0;
    std::string filename;
    uint64_t insertion_counter = 0;
    uint64_t extraction_counter = 0;
    short read_buffer[BUF_CAPACITY];
    short write_buffer[BUF_CAPACITY];
    int read_buffer_pos = 0;
    int read_buffer_fill = 0;
    int write_buffer_pos = 0;
    uint64_t total_shorts_written = 0; // A debug variable.

    explicit file_based_queue(std::string fname) {
        filename = fname;
        appending = fopen(fname.c_str(), "wb");
        reading = fopen(fname.c_str(), "rb");
        setvbuf(reading, NULL, _IONBF, 0);
        assert(reading != nullptr);
        assert(appending != nullptr);

    }

    void fill_read_buffer() {
        size_t succ_read_shorts = 0;
        succ_read_shorts = fread(read_buffer, sizeof(short), BUF_CAPACITY, reading);
        if (succ_read_shorts % (factorial[LISTSIZE]+1) != 0) {
            fprintf(stderr, "Successfully read %zu shorts, but expected to read a multiple of %lu.\n", succ_read_shorts,
                    factorial[LISTSIZE]+1);
            if (feof(reading)) {
                fprintf(stderr, "Reason: end of file.\n");
            } else if (ferror(reading)) {
                fprintf(stderr, "Reason: stream error %d.\n", ferror(reading));
            } else {
                fprintf(stderr, "Reason: not end of file and not a stream error. Hmm...\n");
            }
            abort();
        } else {
            // fprintf(stderr, "Successfully read %zu shorts, 0 mod %lu.\n", succ_read_shorts, factorial[LISTSIZE]+1);
        }
        read_buffer_fill = (int) succ_read_shorts;
        read_buffer_pos = 0;
    }

    void flush_write_buffer() {
        size_t written = fwrite(write_buffer, sizeof(short), write_buffer_pos, appending);
        if (written != write_buffer_pos) {
            PRINT_AND_ABORT("Buffered write failed: not all shorts got written.\n");
        } else {
            // fprintf(stderr, "Written %d (%lu mod %lu) shorts from the buffer.\n", write_buffer_pos,
            //        write_buffer_pos % (factorial[LISTSIZE]+1), factorial[LISTSIZE]+1);
        }
        write_buffer_pos = 0;
        fflush(appending);
    }

    uint64_t size() const {
        return insertion_counter - extraction_counter;
    }

    [[nodiscard]] bool empty() const {
        return insertion_counter <= extraction_counter;
    }

    void push(const workfunction<LISTSIZE>* wf) {
        insertion_counter++;
        assert(write_buffer_pos <= BUF_CAPACITY);
        if (write_buffer_pos == BUF_CAPACITY) {
            flush_write_buffer();
        }
        assert(write_buffer_pos <= BUF_CAPACITY-(factorial[LISTSIZE]+1));
        wf->buffer_serialized_write(write_buffer, write_buffer_pos);
        write_buffer_pos += factorial[LISTSIZE]+1;
        // wf->serialized_write(appending);
    }

    workfunction<LISTSIZE>* pop_nonblocking() {
        if (insertion_counter <= extraction_counter) {
            return nullptr;
        }

        assert(read_buffer_pos <= read_buffer_fill);

        if (read_buffer_pos == read_buffer_fill) {
            if (size() > BUF_WORKFUNCTIONS) {
                fill_read_buffer();
            } else {
                flush_write_buffer();
                fill_read_buffer();
            }
        }

        extraction_counter++;
        auto ret = workfunction<LISTSIZE>::buffer_serialized_read(read_buffer, read_buffer_pos);
        read_buffer_pos += factorial[LISTSIZE]+1;
        return ret;
    }
};