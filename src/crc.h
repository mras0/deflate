#ifndef DEFLATE_CRC_H
#define DEFLATE_CRC_H

#include <stdint.h>

namespace deflate {

uint32_t update_crc32(uint32_t crc, const uint8_t* beg, const uint8_t* end);

} // namespace deflate

#endif
