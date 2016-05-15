#include <string>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>

#include "deflate.h"
#include "crc.h"

using namespace deflate;

#include <fstream>
std::vector<uint8_t> read_file(const std::string& filename)
{
    std::ifstream in(filename, std::ios_base::binary);
    if (!in) throw std::runtime_error(filename + " not found");
    in.seekg(0, std::ios_base::end);
    const auto file_size = static_cast<int>(in.tellg());
    in.seekg(0, std::ios_base::beg);
    std::vector<uint8_t> input(file_size);
    in.read(reinterpret_cast<char*>(&input[0]), input.size());
    return input;
}

std::vector<uint8_t> gunzip(const std::string& filename)
{
    const auto input     = read_file(filename);
    const int  file_size = static_cast<int>(input.size());
    if (file_size < 18) throw std::runtime_error(filename + " is too small to be a gzip file");


    auto invalid = [&filename] () { throw std::runtime_error(filename + " is not a valid gzip file"); };

    // +---+---+---+---+---+---+---+---+---+---+
    // |ID1|ID2|CM |FLG|     MTIME      |XFL|OS|
    // +---+---+---+---+---+---+---+---+---+---+
    if (input[0] != 0x1f || // ID1
        input[1] != 0x8b || // ID2
        input[2] != 8) {    // CM
        invalid();
    }
    enum { FTEXT=1, FHCRC=2, FEXTRA=4, FNAME=8, FCOMMENT=16 };
    const auto flg = input[3];
    int pos = 10;
    if (flg & FEXTRA) {
        const int xlen = input[pos] | (input[pos+1] << 8);
        // TODO: Check whether it's valid to skip this many bytes
        pos += 2 + xlen;
    }
    auto get_zero_terminated_string = [&] () -> std::string {
        std::string res;
        while (pos < file_size) {
            const char c = input[pos++];
            if (!c) return res;
            res += c;
        }
        invalid();
        return {};
    };
    if (flg & FNAME) {
        get_zero_terminated_string(); // filename
    }
    if (flg & FCOMMENT) {
        get_zero_terminated_string(); // comment
    }
    if (flg & FHCRC) {
        pos += 2; // skip 16-bit header CRC
    }
    assert(file_size >= pos + 8);

    int f = file_size - 8;
    const auto crc32 = static_cast<uint32_t>(input[f + 0]) |
        (static_cast<uint32_t>(input[f + 1]) << 8) |
        (static_cast<uint32_t>(input[f + 2]) << 16) |
        (static_cast<uint32_t>(input[f + 3]) << 24);

    const auto isize = static_cast<uint32_t>(input[f + 4]) |
        (static_cast<uint32_t>(input[f + 5]) << 8) |
        (static_cast<uint32_t>(input[f + 6]) << 16) |
        (static_cast<uint32_t>(input[f + 7]) << 24);

    bit_stream bs(input.data() + pos, input.data() + file_size - 8);
    const auto output = deflate::deflate(bs);
    if (output.size() != isize) {
        invalid();
    }

    if (update_crc32(0, output.data(), output.data() + output.size()) != crc32) {
        invalid();
    }

    return output;
}

template<typename F>
void time_it(F f)
{
    constexpr int num_timings = 20;
    double timings[num_timings];
    double sum = 0.0;
    for (int i = 0; i < num_timings; ++i) {
        const auto start = std::chrono::high_resolution_clock::now();
        f();
        const auto end = std::chrono::high_resolution_clock::now();
        timings[i] = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start).count();
        std::cout << timings[i] << "\n";
        sum += timings[i];
    }
    std::sort(timings, timings + num_timings);
    std::cout << "Min/Avg/Mean/Max: ";
    std::cout << timings[0] << " / ";
    std::cout << sum/num_timings << " / ";
    std::cout << timings[num_timings/2] << " / ";
    std::cout << timings[num_timings-1] << std::endl;
}

void timing()
{
    // Tests performned on bunny.tar.gz (4.894.286 B) and doing 2 "warp up" runs before sampling
    //
    // Before optimizations                             Min/Avg/Mean/Max: 245.635 / 262.008 / 249.782 / 347.802
    // Use tables:                                      Min/Avg/Mean/Max: 164.282 / 170.603 / 168.520 / 204.303
    // Remember tables, resize before main deflate loop Min/Avg/Mean/Max: 146.035 / 151.077 / 149.078 / 181.598
    // Rrwrite copy_match to use pointers               Min/Avg/Mean/Max: 144.599 / 149.110 / 146.850 / 178.669
    // Use memcpy in copy_match when possible           Min/Avg/Mean/Max: 140.990 / 145.125 / 143.224 / 166.766
    auto data = read_file("../bunny.tar.gz");
    time_it([&data] {
        bit_stream bs{data.data() + 20, data.data() + data.size() - 8};
        deflate::deflate(bs);
    });
}

int main()
{
    try {
        gunzip("../CMakeLists.txt.gz");
        gunzip("../main.cpp.gz");
        timing();
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
