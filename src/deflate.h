#ifndef DEFLATE_DEFLATE_H
#define DEFLATE_DEFLATE_H

#include "bit_stream.h"
#include <vector>

namespace deflate {

std::vector<uint8_t> deflate(bit_stream& bs);

} // namespace deflate

#endif
