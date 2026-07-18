#pragma once

#include <algorithm>
#include <cmath>

namespace coalescent {
namespace gendyn {

constexpr int kMaxBreakpoints = 64;
constexpr float kDurationErrorMin = -4.f;
constexpr float kDurationErrorMax = 4.f;
constexpr float kLockCorrectionMin = 0.25f;
constexpr float kLockCorrectionMax = 4.f;
constexpr float kLockStepMin = 0.8f;
constexpr float kLockStepMax = 1.25f;

struct DurationBarriers {
    float center;
    float minimum;
    float maximum;
};

struct LockControllerState {
    float correction = 1.f;
    float durationError = 0.f;
    float target = 0.f;
    float normalization = 1.f;
};

// Fold a finite value into [minimum, maximum]. A collapsed interval is a
// deliberate, useful state: DUR WID=0 pins every duration to its center.
inline float reflect(float value, float minimum, float maximum) {
    if (maximum <= minimum)
        return minimum;
    const float range = maximum - minimum;
    float folded = std::fmod(value - minimum, 2.f * range);
    if (folded < 0.f)
        folded += 2.f * range;
    return minimum + (folded > range ? 2.f * range - folded : folded);
}

inline DurationBarriers durationBarriers(float width, float sampleRate,
                                          float centerFrequency, int breakpoints) {
    const float center = sampleRate / (centerFrequency * static_cast<float>(breakpoints));
    const float halfWidth = width * center;
    const float minimum = std::max(1.f, center - halfWidth);
    const float maximum = std::max(minimum, center + halfWidth);
    return {center, minimum, maximum};
}

// One point of the persistent reflected walk. Supplying the random draw keeps
// RNG ownership in the Rack wrapper while making the state transition testable.
inline void advanceWalkPoint(float& value, float& persistentStep, float draw,
                             float gain, float minimum, float maximum) {
    persistentStep = reflect(persistentStep + draw, -1.f, 1.f);
    value = reflect(value + gain * persistentStep, minimum, maximum);
}

inline float durationWalkGain(float scale, const DurationBarriers& barriers) {
    const float halfWidth = (barriers.maximum - barriers.minimum) * 0.5f;
    return scale * halfWidth;
}

// Convert one floating-point segment duration to samples while carrying its
// rounding residual into the next segment. This is the production playback
// quantizer, including the one-sample physical floor and bounded residual.
inline int quantizeDuration(float duration, float normalization, float& error) {
    const float scaled = duration * normalization + error;
    const int samples = std::max(1, static_cast<int>(scaled + 0.5f));
    error = std::fmax(std::fmin(scaled - static_cast<float>(samples),
                                kDurationErrorMax),
                      kDurationErrorMin);
    return samples;
}

inline float lockNormalization(float durationSum, float sampleRate,
                               float centerFrequency, float correction) {
    return (sampleRate / centerFrequency) / durationSum * correction;
}

inline float targetPeriod(bool lockEnabled, float sampleRate, float centerFrequency) {
    return lockEnabled ? sampleRate / centerFrequency : 0.f;
}

// Feedback is intentionally disabled at and below the N-sample floor. There the
// requested period is unreachable and integrating would only store stale error.
inline float applyLockFeedback(float correction, float target, float realized,
                               int breakpoints) {
    if (target > static_cast<float>(breakpoints) && realized > 0.f) {
        const float step = std::fmax(std::fmin(target / realized, kLockStepMax),
                                     kLockStepMin);
        correction = std::fmax(std::fmin(correction * step, kLockCorrectionMax),
                               kLockCorrectionMin);
    }
    return correction;
}

inline bool lockControllerNeedsReset(float oldTarget, float newTarget,
                                     int breakpoints) {
    return newTarget != oldTarget || newTarget <= static_cast<float>(breakpoints);
}

// Install the feed-forward state for a cycle without consuming or changing the
// duration-error remainder. This is used for fresh and restored oscillators.
inline void initializeLockController(LockControllerState& state,
                                     bool lockEnabled, float sampleRate,
                                     float centerFrequency, float durationSum) {
    state.normalization = 1.f;
    if (lockEnabled && durationSum > 0.f) {
        state.normalization = lockNormalization(
            durationSum, sampleRate, centerFrequency, state.correction);
    } else {
        state.correction = 1.f;
    }
    state.target = targetPeriod(lockEnabled, sampleRate, centerFrequency);
}

// Complete the old cycle and install the next one in production order. The
// callback advances the breakpoint walk and returns its new duration sum; RNG
// ownership and draw order therefore remain in the Rack wrapper.
template <typename PrepareNextCycle>
inline void completeCycleAndRetarget(LockControllerState& state,
                                     float realizedPeriod, int breakpoints,
                                     bool nextLockEnabled, float sampleRate,
                                     float nextCenterFrequency,
                                     PrepareNextCycle&& prepareNextCycle) {
    state.correction = applyLockFeedback(
        state.correction, state.target, realizedPeriod, breakpoints);

    const float nextDurationSum = prepareNextCycle();
    const float nextTarget = targetPeriod(
        nextLockEnabled, sampleRate, nextCenterFrequency);
    if (lockControllerNeedsReset(state.target, nextTarget, breakpoints)) {
        state.correction = 1.f;
        state.durationError = 0.f;
    }

    initializeLockController(state, nextLockEnabled, sampleRate,
                             nextCenterFrequency, nextDurationSum);
}

} // namespace gendyn
} // namespace coalescent
