#include <string>
#include <filesystem>
#include <cstdarg>
#include <cassert>
#include <array>
#include <cinttypes>

// A fprintf(stderr) version of the #error macro.
void PRINT_AND_ABORT(const char *format, ...) {
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stderr, format, argptr);
    va_end(argptr);

    abort();
}

std::array<short, 3>* old_format = nullptr;
std::array<int8_t, 3>* new_format = nullptr;
uint64_t advsize = 0;

void init_arrays() {
    old_format = new std::array<short, 3>[advsize];
    new_format = new std::array<int8_t, 3>[advsize];
    for (int i = 0; i < advsize; i++) {
      for (int j = 0; j < 3; j++) {
        old_format[i][j] = 0;
        new_format[i][j] = 0;
      }
    }
}

void deserialize_last_three_short(const std::string& last_three_filename) {

    FILE* binary_file = fopen(last_three_filename.c_str(), "rb");
    size_t read = 0;
    uint64_t advsize_check = 0;
    read = fread(&advsize_check, sizeof(uint64_t), 1, binary_file);
    if (read != 1) {
        PRINT_AND_ABORT("ADVSIZE was not read correctly.");
    }

    fprintf(stderr, "advsize read as %" PRIu64 "\n", advsize_check);
    advsize = advsize_check;

    init_arrays();

    read = fread(old_format, sizeof(std::array<short, 3>), advsize, binary_file);
    if (read != advsize) {
        PRINT_AND_ABORT("The last three choices array was not read correctly.");
    }

    fclose(binary_file);

    fprintf(stderr, "0: %hd, %hd, %hd.\n",
        old_format[0][0], old_format[0][1], old_format[0][2]);
    fprintf(stderr, "1: %hd, %hd, %hd.\n",
        old_format[1][0], old_format[1][1], old_format[2][2]);
    fprintf(stderr, "2: %hd, %hd, %hd.\n",
        old_format[2][0], old_format[2][1], old_format[2][2]);
}

void sync_arrays() {
    for (int i = 0; i < advsize; i++) {
        for (int j = 0; j < 3; j++) {
          new_format[i][j] = static_cast<int8_t>(old_format[i][j]);
        }
    }
}

void serialize_last_three_int8t(const std::string& last_three_filename) {

        FILE* binary_file = fopen(last_three_filename.c_str(), "wb");
        size_t written = 0;
        written = fwrite(&advsize, sizeof(uint64_t), 1, binary_file);
        if (written != 1) {
            PRINT_AND_ABORT("ADVSIZE was not written correctly.");
        }
        written = fwrite(new_format, sizeof(std::array<int8_t, 3>), advsize, binary_file);
        if (written != advsize) {
            PRINT_AND_ABORT("The last three choices array was not written correctly.");
        }

        fclose(binary_file);

        fprintf(stderr, "0: %" PRIi8 ", %" PRIi8 ", %" PRIi8 ".\n",
            new_format[0][0], new_format[0][1], new_format[0][2]);
        fprintf(stderr, "1: %" PRIi8 ", %" PRIi8 ", %" PRIi8 ".\n",
            new_format[1][0], new_format[1][1], new_format[2][2]);
        fprintf(stderr, "2: %" PRIi8 ", %" PRIi8 ", %" PRIi8 ".\n",
            new_format[2][0], new_format[2][1], new_format[2][2]);
    }

int main(void)
{
    std::string last_three_filename = std::string("last-three-maximizers-5-shorts.bin");
    std::string last_three_filename_int8t = std::string("last-three-maximizers-5-int8t.bin");

    deserialize_last_three_short(last_three_filename);
    sync_arrays();
    serialize_last_three_int8t(last_three_filename_int8t);
    return 0;
}