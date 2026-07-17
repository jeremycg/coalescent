#pragma once

#include "../dsp/pcg32.hpp"

namespace coalescent {
namespace lineages {

// Preserve the Lineages-facing name while sharing the exact PCG engine with
// Islands. Generated sequences and serialized state remain unchanged.
using coalescent::Pcg32;

} // namespace lineages
} // namespace coalescent
