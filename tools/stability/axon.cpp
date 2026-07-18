// Regression and calibration checks for Axon's FitzHugh-Nagumo model.
//
// The transition, constants, scheduling, and safety repair below are the exact
// SDK-free production core used by src/neuron/Axon.cpp. This test deliberately
// keeps its measurements independent: threshold crossing, period estimation,
// pitch error, and finite/bounds assertions are external oracles around that
// shared state transition.
#include "../../src/dsp/neuron_models.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

using Core = coalescent::neuron::AxonCore;

static inline bool upwardSpikeCrossing(float previous, float current) {
    return previous < Core::SPIKE_THRESH && current >= Core::SPIKE_THRESH;
}

static inline bool rawStateIsSafe(const float (&state)[Core::STATE_COUNT]) {
    for (int i = 0; i < Core::STATE_COUNT; ++i) {
        if (!std::isfinite(state[i]) || std::fabs(state[i]) >= Core::RUNAWAY_MAX)
            return false;
    }
    return true;
}

static inline bool rawStateIsWithinClamp(const float (&state)[Core::STATE_COUNT]) {
    for (int i = 0; i < Core::STATE_COUNT; ++i) {
        if (std::fabs(state[i]) > Core::STATE_MAX + 1e-3f)
            return false;
    }
    return true;
}

static inline void resetState(float (&state)[Core::STATE_COUNT]) {
    state[0] = Core::REST_V;
    state[1] = Core::REST_W;
}

static inline bool advanceScheduled(float (&state)[Core::STATE_COUNT],
                                    const coalescent::neuron::ScalarSchedule& schedule,
                                    float current, float eps, float shape) {
    Core::advanceObservation(state, schedule.h, schedule.substeps, current, eps, shape);
    const bool rawSafe = rawStateIsSafe(state);
    const bool retained = Core::repair(state);
    return rawSafe && retained;
}

// Fine model-time reference. The independent crossing/timing oracle sees state
// only at the same accepted observation boundary as the Rack wrapper.
static double measurePeriodTau(float current, float eps, float shape) {
    float state[Core::STATE_COUNT];
    resetState(state);
    const float dstep = 1e-4f;
    const coalescent::neuron::ScalarSchedule schedule =
        coalescent::neuron::scalarScheduleFromSubTau<Core>(dstep);

    double tau = 0.0;
    for (long i = 0; i < (long) (400.0 / dstep); ++i) {
        advanceScheduled(state, schedule, current, eps, shape);
        tau += dstep;
    }

    float previous = state[0];
    double firstCross = -1.0;
    double lastCross = -1.0;
    int crossings = 0;
    const long maxIterations = (long) (2000.0 / dstep);
    for (long i = 0; i < maxIterations && crossings < 21; ++i) {
        advanceScheduled(state, schedule, current, eps, shape);
        tau += dstep;
        if (upwardSpikeCrossing(previous, state[0])) {
            if (firstCross < 0.0)
                firstCross = tau;
            lastCross = tau;
            ++crossings;
        }
        previous = state[0];
    }
    return crossings >= 2 ? (lastCross - firstCross) / (crossings - 1) : -1.0;
}

int main() {
    int failures = 0;

    // Production constants must still agree with an independently measured
    // default limit-cycle period; the model is open-loop, but calibration itself
    // should not silently drift when the shared core changes.
    const double periodTau = measurePeriodTau(0.6f, 0.08f, 0.7f);
    if (periodTau <= 0.0) {
        std::printf("FAIL: default parameters do not oscillate (T_dim=%.4f)\n", periodTau);
        ++failures;
    }
    else {
        const double calibrationCents = 1200.0 * std::log2(Core::RATE_CAL / periodTau);
        std::printf("CALIBRATION: measured T_dim=%.6f, production RATE_CAL=%.6f (%+.3f cents)\n",
                    periodTau, Core::RATE_CAL, calibrationCents);
        if (std::fabs(calibrationCents) > 1.0) {
            std::printf("  FAIL: production calibration differs from the fine reference by >1 cent\n");
            ++failures;
        }
    }

    // Contract checks around the shared schedule and sanitizers use literal
    // expected results rather than recomputing the production formula locally.
    {
        const coalescent::neuron::ScalarSchedule capped =
            coalescent::neuron::scalarSchedule<Core>(1e30f, 1.f);
        const float cap = Core::HSUB_MAX * Core::MAX_SUB;
        const bool scheduleOk = capped.subTau == cap
            && capped.substeps == Core::MAX_SUB
            && std::fabs(capped.h - Core::HSUB_MAX) < 1e-7f;
        const bool pitchGuardOk = coalescent::neuron::sanitizePitchExponent(
                std::numeric_limits<float>::quiet_NaN()) == 0.f
            && coalescent::neuron::sanitizePitchExponent(
                std::numeric_limits<float>::infinity()) == 30.f
            && coalescent::neuron::sanitizePitchExponent(
                -std::numeric_limits<float>::infinity()) == -30.f;
        std::printf("SCHEDULE cap/substeps and hostile pitch guard: %s\n",
                    scheduleOk && pitchGuardOk ? "PASS" : "FAIL");
        if (!scheduleOk || !pitchGuardOk)
            ++failures;
    }

    // The bounded-CPU policy is user-visible: at 44.1 kHz the simulation-speed
    // ceiling moves up by one octave for each doubling of oversampling. Keep
    // literal brackets around the documented limits so a scheduler change
    // cannot silently move the playable range.
    {
        struct CeilingBracket {
            int oversample;
            float belowOctaves;
            float aboveOctaves;
        };
        const CeilingBracket brackets[] = {
            {1, 3.82f, 3.84f},
            {2, 4.82f, 4.84f},
            {4, 5.82f, 5.84f},
            {8, 6.82f, 6.84f},
        };
        bool ceilingsOk = true;
        const float cap = Core::HSUB_MAX * Core::MAX_SUB;
        for (const CeilingBracket& bracket : brackets) {
            const coalescent::neuron::ScalarSchedule below =
                coalescent::neuron::scalarSchedule<Core>(
                    coalescent::neuron::PITCH_REFERENCE_HZ
                        * std::exp2(bracket.belowOctaves),
                    44100.f, bracket.oversample);
            const coalescent::neuron::ScalarSchedule above =
                coalescent::neuron::scalarSchedule<Core>(
                    coalescent::neuron::PITCH_REFERENCE_HZ
                        * std::exp2(bracket.aboveOctaves),
                    44100.f, bracket.oversample);
            ceilingsOk = ceilingsOk && below.subTau < cap && above.subTau == cap;
        }
        std::printf("SCHEDULE documented 44.1 kHz ceilings: %s\n",
                    ceilingsOk ? "PASS" : "FAIL");
        if (!ceilingsOk)
            ++failures;
    }

    // Exercise the exact shared repair contract separately from the normal-run
    // invariant below: huge/NaN states reset atomically, finite overshoot clamps.
    {
        float runaway[Core::STATE_COUNT] = {
            std::numeric_limits<float>::quiet_NaN(), 0.f
        };
        const bool runawayRetained = Core::repair(runaway);
        float finite[Core::STATE_COUNT] = {15.f, -15.f};
        const bool finiteRetained = Core::repair(finite);
        const bool repairOk = !runawayRetained
            && runaway[0] == Core::REST_V && runaway[1] == Core::REST_W
            && finiteRetained && finite[0] == Core::STATE_MAX
            && finite[1] == -Core::STATE_MAX;
        std::printf("SAFETY reset/clamp contract: %s\n", repairOk ? "PASS" : "FAIL");
        if (!repairOk)
            ++failures;
    }

    // CURRENT CV is neutral lane-by-lane when a cable supplies NaN/Inf. The
    // generic shared policy is instantiated here with scalar predicates; Rack
    // supplies SIMD predicates/selectors around the same entry point.
    {
        const auto finite = [](float value) { return std::isfinite(value); };
        const auto select = [](bool condition, float yes, float no) {
            return condition ? yes : no;
        };
        const bool currentGuardOk = coalescent::finiteOr(2.5f, 0.f, finite, select) == 2.5f
            && coalescent::finiteOr(
                std::numeric_limits<float>::quiet_NaN(), 0.f, finite, select) == 0.f
            && coalescent::finiteOr(
                std::numeric_limits<float>::infinity(), 0.f, finite, select) == 0.f
            && coalescent::finiteOr(
                -std::numeric_limits<float>::infinity(), 0.f, finite, select) == 0.f;
        std::printf("CURRENT CV non-finite neutralization: %s\n",
                    currentGuardOk ? "PASS" : "FAIL");
        if (!currentGuardOk)
            ++failures;
    }

    const float sampleRate = 48000.f;
    int sweepCount = 0;
    bool rawClampUsed = false;
    for (float current = -0.2f; current <= 1.6001f; current += 0.2f) {
        for (float eps = 0.01f; eps <= 0.3001f; eps += 0.04f) {
            for (float shape = 0.4f; shape <= 1.0001f; shape += 0.15f) {
                for (float oct = -4.f; oct <= 4.001f; oct += 4.f) {
                    float state[Core::STATE_COUNT];
                    resetState(state);
                    const float pitchHz = coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
                    const coalescent::neuron::ScalarSchedule schedule =
                        coalescent::neuron::scalarSchedule<Core>(pitchHz, sampleRate);
                    const long samples = (long) (sampleRate * 0.5f);
                    for (long i = 0; i < samples; ++i) {
                        Core::advanceObservation(
                            state, schedule.h, schedule.substeps, current, eps, shape);
                        const bool rawSafe = rawStateIsSafe(state);
                        const bool withinClamp = rawStateIsWithinClamp(state);
                        rawClampUsed = rawClampUsed || !withinClamp;
                        const bool retained = Core::repair(state);
                        if (!rawSafe || !retained || !rawStateIsWithinClamp(state)) {
                            std::printf("FAIL: state escaped at I=%.2f eps=%.3f a=%.2f oct=%.0f: v=%g w=%g\n",
                                        current, eps, shape, oct, state[0], state[1]);
                            ++failures;
                            break;
                        }
                    }
                    ++sweepCount;
                }
            }
        }
    }
    std::printf("STABILITY: %d parameter/pitch combinations finite and bounded; raw clamp %s\n",
                sweepCount, rawClampUsed ? "was needed" : "was not needed");

    // Independent event timing over the production schedule. This intentionally
    // measures the emergent oscillator, not an algebraic restatement of RATE_CAL.
    auto measureAudioHz = [&](float oct, float current, float eps, float shape) -> double {
        float state[Core::STATE_COUNT];
        resetState(state);
        const float pitchHz = coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
        const coalescent::neuron::ScalarSchedule schedule =
            coalescent::neuron::scalarSchedule<Core>(pitchHz, sampleRate);
        for (long i = 0; i < (long) (sampleRate * 0.3f); ++i)
            advanceScheduled(state, schedule, current, eps, shape);

        float previous = state[0];
        double firstCross = -1.0;
        double lastCross = -1.0;
        int crossings = 0;
        const long maxSamples = (long) (sampleRate * 2.f);
        for (long sample = 0; sample < maxSamples && crossings < 21; ++sample) {
            advanceScheduled(state, schedule, current, eps, shape);
            if (upwardSpikeCrossing(previous, state[0])) {
                if (firstCross < 0.0)
                    firstCross = sample;
                lastCross = sample;
                ++crossings;
            }
            previous = state[0];
        }
        if (crossings < 2)
            return -1.0;
        return sampleRate / ((lastCross - firstCross) / (crossings - 1));
    };

    for (float oct = -2.f; oct <= 2.001f; oct += 1.f) {
        const double wanted = coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
        const double measured = measureAudioHz(oct, 0.6f, 0.08f, 0.7f);
        const double cents = measured > 0.0 ? 1200.0 * std::log2(measured / wanted) : 0.0;
        std::printf("PITCH: oct=%+.0f want=%7.2f Hz got=%7.2f Hz (%+6.1f cents)\n",
                    oct, wanted, measured, cents);
        if (measured < 0.0 || std::fabs(cents) > 50.0) {
            std::printf("  FAIL: pitch off by more than 50 cents\n");
            ++failures;
        }
    }

    if (failures) {
        std::printf("\n%d CHECK(S) FAILED\n", failures);
        return 1;
    }
    std::printf("\nAll checks passed.\n");
    return 0;
}
