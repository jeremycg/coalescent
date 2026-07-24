#include "../../src/dsp/finite.hpp"
#include "../../src/dsp/neuron_models.hpp"

#include <cstdio>
#include <cstdint>
#include <cstring>

namespace {

float floatFromBits(std::uint32_t bits) {
    volatile std::uint32_t source = bits;
    const std::uint32_t loaded = source;
    float value = 0.f;
    std::memcpy(&value, &loaded, sizeof(value));
    return value;
}

double doubleFromBits(std::uint64_t bits) {
    volatile std::uint64_t source = bits;
    const std::uint64_t loaded = source;
    double value = 0.0;
    std::memcpy(&value, &loaded, sizeof(value));
    return value;
}

void expect(bool condition, const char* contract, int& failures) {
    if (!condition) {
        std::printf("FAIL: %s\n", contract);
        ++failures;
    }
}

} // namespace

int main() {
    int failures = 0;

    const float positiveZero = floatFromBits(UINT32_C(0x00000000));
    const float negativeZero = floatFromBits(UINT32_C(0x80000000));
    const float subnormal = floatFromBits(UINT32_C(0x00000001));
    const float maximum = floatFromBits(UINT32_C(0x7f7fffff));
    const float positiveInfinity = floatFromBits(UINT32_C(0x7f800000));
    const float negativeInfinity = floatFromBits(UINT32_C(0xff800000));
    const float positiveNaN = floatFromBits(UINT32_C(0x7fc12345));
    const float negativeNaN = floatFromBits(UINT32_C(0xffc12345));

    expect(coalescent::isFinite(positiveZero)
               && coalescent::isFinite(negativeZero)
               && coalescent::isFinite(subnormal)
               && coalescent::isFinite(maximum)
               && !coalescent::isFinite(positiveInfinity)
               && !coalescent::isFinite(negativeInfinity)
               && !coalescent::isFinite(positiveNaN)
               && !coalescent::isFinite(negativeNaN),
           "float finite classification", failures);
    expect(!coalescent::isNaN(positiveZero)
               && !coalescent::isNaN(subnormal)
               && !coalescent::isNaN(maximum)
               && !coalescent::isNaN(positiveInfinity)
               && !coalescent::isNaN(negativeInfinity)
               && coalescent::isNaN(positiveNaN)
               && coalescent::isNaN(negativeNaN),
           "float NaN classification", failures);

    const double doubleMaximum =
        doubleFromBits(UINT64_C(0x7fefffffffffffff));
    const double doubleSubnormal =
        doubleFromBits(UINT64_C(0x0000000000000001));
    const double doubleInfinity =
        doubleFromBits(UINT64_C(0x7ff0000000000000));
    const double doubleNegativeInfinity =
        doubleFromBits(UINT64_C(0xfff0000000000000));
    const double doubleNaN =
        doubleFromBits(UINT64_C(0x7ff8123456789abc));
    const double doubleNegativeNaN =
        doubleFromBits(UINT64_C(0xfff8123456789abc));

    expect(coalescent::isFinite(doubleMaximum)
               && coalescent::isFinite(doubleSubnormal)
               && !coalescent::isFinite(doubleInfinity)
               && !coalescent::isFinite(doubleNegativeInfinity)
               && !coalescent::isFinite(doubleNaN)
               && !coalescent::isFinite(doubleNegativeNaN),
           "double finite classification", failures);
    expect(!coalescent::isNaN(doubleMaximum)
               && !coalescent::isNaN(doubleInfinity)
               && !coalescent::isNaN(doubleNegativeInfinity)
               && coalescent::isNaN(doubleNaN)
               && coalescent::isNaN(doubleNegativeNaN),
           "double NaN classification", failures);

    expect(coalescent::finiteOr(2.5f, -7.f) == 2.5f
               && coalescent::finiteOr(positiveNaN, -7.f) == -7.f
               && coalescent::finiteOr(positiveInfinity, -7.f) == -7.f
               && coalescent::finiteOr(negativeInfinity, -7.f) == -7.f,
           "scalar finite fallback", failures);

    expect(coalescent::neuron::sanitizePitchExponent(2.5f) == 2.5f
               && coalescent::neuron::sanitizePitchExponent(maximum) == 30.f
               && coalescent::neuron::sanitizePitchExponent(-maximum) == -30.f
               && coalescent::neuron::sanitizePitchExponent(positiveInfinity) == 30.f
               && coalescent::neuron::sanitizePitchExponent(negativeInfinity) == -30.f
               && coalescent::neuron::sanitizePitchExponent(positiveNaN) == 0.f
               && coalescent::neuron::sanitizePitchExponent(negativeNaN) == 0.f,
           "scalar neuron pitch sanitization", failures);

    float axonState[coalescent::neuron::AxonCore::STATE_COUNT] = {
        positiveNaN, 0.f
    };
    const bool axonRetained =
        coalescent::neuron::AxonCore::repair(axonState);
    expect(!axonRetained
               && axonState[0] == coalescent::neuron::AxonCore::REST_V
               && axonState[1] == coalescent::neuron::AxonCore::REST_W,
           "Axon scalar state repair", failures);

    float somaState[coalescent::neuron::SomaCore::STATE_COUNT] = {
        0.f, negativeInfinity, 0.f
    };
    const bool somaRetained =
        coalescent::neuron::SomaCore::repair(somaState);
    expect(!somaRetained
               && somaState[0] == coalescent::neuron::SomaCore::REST_X
               && somaState[1] == coalescent::neuron::SomaCore::REST_Y
               && somaState[2] == coalescent::neuron::SomaCore::REST_Z,
           "Soma scalar state repair", failures);

    if (failures)
        return 1;

    std::printf("finite classification and neuron sanitizers: PASS\n");
    return 0;
}
