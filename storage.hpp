#pragma once

#include <cstdio>
#include "common.hpp"

void write_distance_array(cost_t *array, unsigned long int size, FILE *write_file = nullptr) {
    bool predefined_file_mode = (write_file == nullptr);
    if (predefined_file_mode) {
        char filename[256] = {0};
        sprintf(filename, "./distance-file-%d.bin", LISTSIZE);
        write_file = fopen(filename, "wb");
    }

    fwrite(&size, sizeof(unsigned long int), 1, write_file);
    fwrite(array, sizeof(cost_t), size, write_file);

    if (predefined_file_mode) {
        fclose(write_file);
    }
}

std::pair<unsigned long int, cost_t*> read_distance_array(FILE *read_file = nullptr) {
    bool predefined_file_mode = (read_file == nullptr);
    if (predefined_file_mode) {
        char filename[256] = {0};
        sprintf(filename, "./distance-file-%d.bin", LISTSIZE);
        read_file = fopen(filename, "rb");
    }

    unsigned long int distance_array_length = 0;
    size_t length_check = fread(&distance_array_length, sizeof(unsigned long int), 1, read_file);
    if (length_check != 1) {
        fprintf(stderr, "Failed to read the length of the file.\n");
        exit(-1);
    }
    cost_t *distances = new cost_t[distance_array_length];
    size_t elements_read = fread(distances, sizeof(cost_t), distance_array_length, read_file);

    if (elements_read != distance_array_length) {
        fprintf(stderr, "Read fewer elements than advertised.\n");
        exit(-1);
    }
    if (predefined_file_mode) {
        fclose(read_file);
    }

    return {distance_array_length, distances};
}