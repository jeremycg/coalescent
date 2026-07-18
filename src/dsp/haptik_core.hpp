#pragma once

#include "ode_policy.hpp"

#include <algorithm>
#include <cmath>

namespace coalescent {
namespace haptik {

constexpr int kMaxStorageMasses = 256;
constexpr int kMinPlayableMasses = 8;
constexpr int kMaxPlayableMasses = 128;
constexpr int kSlowDivider = 256;

constexpr float kPi = 3.14159265358979323846f;
constexpr float kDampingMaxHz = 250.f;
constexpr float kOutputGain = 1.f;
constexpr float kMotionGain = 4.f;
constexpr float kExternalGain = 0.005f;
constexpr float kBumpFraction = 0.125f;
constexpr float kCvDepth = 0.1f;
constexpr float kDriveKeepalive = 0.05f;
constexpr float kStateMax = 16.f;
constexpr float kCenteringMax = 0.35f;
constexpr float kCouplingMax = 0.9f;
constexpr float kEvolutionMinHz = 0.05f;
constexpr float kEvolutionMaxHz = 30.f;
constexpr float kDenormalGuard = 1e-20f;
constexpr float kLog2E = 1.44269504f;

inline float externalDrive(float voltageSum, float amount) {
    return coalescent::finiteOr(voltageSum, 0.f) * kExternalGain * amount;
}

inline float centeringCoefficient(float evolutionHz, int divider, float sampleRate) {
    const float omega = 2.f * kPi * evolutionHz * divider / sampleRate;
    return std::min(omega * omega, kCenteringMax);
}

// Argument for the Rack exp2 approximation used by the wrapper. Keeping the
// approximation itself in Rack avoids vendoring SDK math into this SDK-free core.
inline float dampingExp2Argument(float damping, int divider, float sampleRate) {
    const float tapered = damping * damping;
    return -tapered * kDampingMaxHz * divider / sampleRate * kLog2E;
}

inline float nextVelocity(float velocity, float left, float center, float right,
                          float coupling, float centering, float dampingMultiplier) {
    const float acceleration = coupling * (left - 2.f * center + right)
                             - centering * center;
    return (velocity + acceleration) * dampingMultiplier + kDenormalGuard;
}

inline void advanceVelocities(const float* displacement, float* velocity, int masses,
                              float coupling, float centering,
                              float dampingMultiplier) {
    velocity[0] = nextVelocity(velocity[0], displacement[masses - 1], displacement[0],
                               displacement[1], coupling, centering, dampingMultiplier);
    for (int i = 1; i < masses - 1; ++i) {
        velocity[i] = nextVelocity(velocity[i], displacement[i - 1], displacement[i],
                                   displacement[i + 1], coupling, centering,
                                   dampingMultiplier);
    }
    velocity[masses - 1] = nextVelocity(velocity[masses - 1], displacement[masses - 2],
                                        displacement[masses - 1], displacement[0],
                                        coupling, centering, dampingMultiplier);
}

inline void applyDriverForce(float& velocity, float drive, float dampingMultiplier) {
    velocity += drive * dampingMultiplier;
}

// The callable min/max arguments let the same operation order serve scalar
// std::fmin/std::fmax and Rack's SIMD equivalents without depending on the SDK.
template <typename T, typename Min, typename Max>
inline T clampStateValue(T value, const T& minimum, const T& maximum,
                         Min&& minimumOf, Max&& maximumOf) {
    return maximumOf(minimumOf(value, maximum), minimum);
}

template <typename T, typename Min, typename Max>
inline void advanceStateValue(T& displacement, T& velocity,
                              const T& minimum, const T& maximum,
                              Min&& minimumOf, Max&& maximumOf) {
    displacement = clampStateValue(displacement + velocity, minimum, maximum,
                                   minimumOf, maximumOf);
    velocity = clampStateValue(velocity, minimum, maximum,
                               minimumOf, maximumOf);
}

inline float clampStateValue(float value) {
    return std::fmax(std::fmin(value, kStateMax), -kStateMax);
}

inline void advanceStateValue(float& displacement, float& velocity) {
    displacement = clampStateValue(displacement + velocity);
    velocity = clampStateValue(velocity);
}

inline void stepScalar(float* displacement, float* velocity, int masses,
                       int driver, float drive, float coupling, float centering,
                       float dampingMultiplier) {
    advanceVelocities(displacement, velocity, masses, coupling, centering,
                      dampingMultiplier);
    applyDriverForce(velocity[driver], drive, dampingMultiplier);
    for (int i = 0; i < masses; ++i)
        advanceStateValue(displacement[i], velocity[i]);
}

inline float interpolate(float previous, float current, float fraction) {
    return previous + fraction * (current - previous);
}

inline void collapseInterpolatedFrame(const float* previous, float* current,
                                      int masses, float fraction) {
    for (int i = 0; i < masses; ++i)
        current[i] = interpolate(previous[i], current[i], fraction);
}

inline void captureInterpolatedFrame(float* previous, float* current,
                                     int masses, float fraction) {
    collapseInterpolatedFrame(previous, current, masses, fraction);
    std::copy(current, current + masses, previous);
}

inline bool shouldStep(bool frozen, bool slow, int dividerCounter) {
    return !frozen && (!slow || dividerCounter == 0);
}

inline int advanceDivider(int dividerCounter, int divider) {
    return (dividerCounter + 1) % divider;
}

inline void advanceScanPhase(float& phase, float pitchHz, float sampleRate) {
    phase += pitchHz / sampleRate;
    if (!std::isfinite(phase))
        phase = 0.f;
    phase -= std::floor(phase);
}

inline float addImpulse(float* displacement, int driver, float amount) {
    displacement[driver] += amount;
    return amount;
}

inline float addHannBump(float* displacement, int masses, int driver, float amount) {
    float added = 0.f;
    const int halfWidth = std::max(1, static_cast<int>(kBumpFraction * masses));
    for (int distance = -halfWidth; distance <= halfWidth; ++distance) {
        const float weight = 0.5f * (1.f + std::cos(kPi * distance / halfWidth));
        const int index = ((driver + distance) % masses + masses) % masses;
        displacement[index] += amount * weight;
        added += amount * weight;
    }
    return added;
}

inline void removeInjectedMean(float* displacement, int masses, float added) {
    if (added == 0.f)
        return;
    const float mean = added / masses;
    for (int i = 0; i < masses; ++i)
        displacement[i] -= mean;
}

} // namespace haptik
} // namespace coalescent
