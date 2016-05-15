#include "crc.h"
#include <array>
#include <utility>

namespace deflate {

constexpr uint32_t crc32_poly = 0xedb88320; // 0x04C11DB7 reversed

constexpr uint32_t crc32_one_bit(uint32_t c)
{
    return c & 1 ? crc32_poly ^ (c >> 1) : c >> 1;
}

constexpr uint32_t crc32_one_byte(uint32_t c)
{
    return crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(crc32_one_bit(c))))))));
}

#define USE_CRC32_TABLE
#ifdef USE_CRC32_TABLE
template<size_t... is>
constexpr std::array<uint32_t, 256> make_crc32_table(std::index_sequence<is...>)
{
    return { crc32_one_byte(is)... };
}

constexpr std::array<uint32_t, 256> crc32_table = make_crc32_table(std::make_index_sequence<256>());
static_assert(crc32_table[0] == 0, "");
static_assert(crc32_table[255] == 0x2d02ef8d, "");
#endif

uint32_t update_crc32(uint32_t crc, const uint8_t* beg, const uint8_t* end)
{
    crc ^= ~0U;
    for (auto p = beg; p != end; ++p) {
#ifdef USE_CRC32_TABLE
        crc = crc32_table[(crc ^ *p) & 0xff] ^ (crc >> 8);
#else
        crc = crc32_one_byte(crc ^ *p);
#endif
    }
    crc ^= ~0U;
    return crc;
}

} // namespace deflate
