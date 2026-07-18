// Regression and calibration checks for Soma's Hindmarsh-Rose model.
//
// The transition, constants, scheduling, parameter sanitization, and safety
// repair are the exact SDK-free production core used by src/neuron/Soma.cpp.
// Period selection, threshold crossings, pitch error, and finite/bounds checks
// remain independent measurement oracles around that shared transition.
#include "../../src/dsp/neuron_models.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using Core = coalescent::neuron::SomaCore;

static constexpr float DEFAULT_CURRENT = 2.f;
static constexpr float DEFAULT_ADAPT = 4.f;

static inline bool upwardSpikeCrossing(float previous, float current) {
    return previous < Core::SPIKE_THRESH && current >= Core::SPIKE_THRESH;
}

static inline void resetState(float (&state)[Core::STATE_COUNT]) {
    state[0] = Core::REST_X;
    state[1] = Core::REST_Y;
    state[2] = Core::REST_Z;
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

static inline bool advanceScheduled(float (&state)[Core::STATE_COUNT],
                                    const coalescent::neuron::ScalarSchedule& schedule,
                                    float current, float r, float adapt) {
    Core::advanceObservation(state, schedule.h, schedule.substeps, current, r, adapt);
    const bool rawSafe = rawStateIsSafe(state);
    const bool retained = Core::repair(state);
    return rawSafe && retained;
}

// Independent musical oracle: retain within-burst intervals and reject long
// inter-burst gaps before taking the median.
static double medianWithinBurst(std::vector<double> intervals) {
    if (intervals.size() < 3)
        return -1.0;
    std::sort(intervals.begin(), intervals.end());
    const double minimum = intervals.front();
    std::vector<double> within;
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        if (intervals[i] < 2.5 * minimum)
            within.push_back(intervals[i]);
    }
    return within.empty() ? -1.0 : within[within.size() / 2];
}

// Fine fixed-step model-time reference using the shared derivative/RK4 kernel,
// with independent event timing and burst classification.
static double measurePeriodTau(float current, float r, float adapt) {
    float state[Core::STATE_COUNT];
    resetState(state);
    const float step = 5e-4f;
    for (long i = 0; i < (long) (800.0 / step); ++i)
        Core::advanceObservation(state, step, 1, current, r, adapt);

    float previous = state[0];
    std::vector<double> intervals;
    double lastCross = -1.0;
    double tau = 0.0;
    for (long i = 0; i < (long) (3000.0 / step); ++i) {
        Core::advanceObservation(state, step, 1, current, r, adapt);
        tau += step;
        if (upwardSpikeCrossing(previous, state[0])) {
            if (lastCross >= 0.0)
                intervals.push_back(tau - lastCross);
            lastCross = tau;
        }
        previous = state[0];
    }
    return medianWithinBurst(intervals);
}

int main() {
    int failures = 0;

    const double spikePeriod =
        measurePeriodTau(DEFAULT_CURRENT, Core::R_DEFAULT, DEFAULT_ADAPT);
    if (spikePeriod <= 0.0) {
        std::printf("FAIL: default voicing does not spike (T=%.3f)\n", spikePeriod);
        ++failures;
    }
    else {
        const double calibrationCents = 1200.0 * std::log2(Core::RATE_CAL / spikePeriod);
        std::printf("CALIBRATION: measured within-burst T=%.6f, production RATE_CAL=%.6f (%+.3f cents)\n",
                    spikePeriod, Core::RATE_CAL, calibrationCents);
        if (std::fabs(calibrationCents) > 1.0) {
            std::printf("  FAIL: production calibration differs from fine reference by >1 cent\n");
            ++failures;
        }
    }

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

    // Literal public pitch-ceiling contract at 44.1 kHz. Just below each
    // documented boundary the schedule remains open-loop; just above it the
    // model-time step saturates at the bounded CPU/stability ceiling.
    {
        const float sampleRate = 44100.f;
        const int oversampling[] = {1, 2, 4, 8};
        const float ceilingOctaves[] = {
            3.28432668f, 4.28432668f, 5.28432668f, 6.28432668f
        };
        const float cap = Core::HSUB_MAX * Core::MAX_SUB;
        bool ceilingsOk = true;
        for (int i = 0; i < 4; ++i) {
            const float belowPitch = coalescent::neuron::PITCH_REFERENCE_HZ
                * std::exp2(ceilingOctaves[i] - 0.001f);
            const float abovePitch = coalescent::neuron::PITCH_REFERENCE_HZ
                * std::exp2(ceilingOctaves[i] + 0.001f);
            const coalescent::neuron::ScalarSchedule below =
                coalescent::neuron::scalarSchedule<Core>(
                    belowPitch, sampleRate, oversampling[i]);
            const coalescent::neuron::ScalarSchedule above =
                coalescent::neuron::scalarSchedule<Core>(
                    abovePitch, sampleRate, oversampling[i]);
            ceilingsOk = ceilingsOk && below.subTau < cap
                && above.subTau == cap && above.substeps == Core::MAX_SUB
                && std::fabs(above.h - Core::HSUB_MAX) < 1e-7f;
        }
        std::printf("SCHEDULE 44.1 kHz ceilings (+3.28/+4.28/+5.28/+6.28): %s\n",
                    ceilingsOk ? "PASS" : "FAIL");
        if (!ceilingsOk)
            ++failures;
    }

    {
        float runaway[Core::STATE_COUNT] = {
            0.f, std::numeric_limits<float>::infinity(), 0.f
        };
        const bool runawayRetained = Core::repair(runaway);
        float finite[Core::STATE_COUNT] = {30.f, -30.f, 0.f};
        const bool finiteRetained = Core::repair(finite);
        const bool repairOk = !runawayRetained
            && runaway[0] == Core::REST_X && runaway[1] == Core::REST_Y
            && runaway[2] == Core::REST_Z && finiteRetained
            && finite[0] == Core::STATE_MAX && finite[1] == -Core::STATE_MAX
            && finite[2] == 0.f;
        std::printf("SAFETY reset/clamp contract: %s\n", repairOk ? "PASS" : "FAIL");
        if (!repairOk)
            ++failures;
    }

    {
        const auto finite = [](float value) { return std::isfinite(value); };
        const auto select = [](bool condition, float yes, float no) {
            return condition ? yes : no;
        };
        const bool currentGuardOk = coalescent::finiteOr(-4.f, 0.f, finite, select) == -4.f
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

    // Exact shared BURST policy, checked against independent expected outcomes.
    {
        const float base = Core::sanitizeBurstBase(std::log2(Core::R_DEFAULT));
        const bool burstGuardOk = Core::sanitizeBurstBase(
                std::numeric_limits<float>::quiet_NaN()) == base
            && Core::sanitizeBurstExponent(
                std::numeric_limits<float>::quiet_NaN(), base) == base
            && Core::sanitizeBurstExponent(
                std::numeric_limits<float>::infinity(), base) == base
            && Core::sanitizeBurstExponent(
                -std::numeric_limits<float>::infinity(), base) == base
            && Core::sanitizeBurstExponent(1e30f, base) == Core::BURST_EXP_MAX;
        std::printf("BURST hostile exponent guard: %s\n", burstGuardOk ? "PASS" : "FAIL");
        if (!burstGuardOk)
            ++failures;
    }

    // At Soma's maximum accepted observation interval, a narrow spike can rise
    // and fall entirely between endpoint observations. The production callback
    // must expose every accepted RK4 state so the Rack Schmitt trigger sees it.
    {
        const coalescent::neuron::ScalarSchedule schedule =
            coalescent::neuron::scalarScheduleFromSubTau<Core>(
                Core::HSUB_MAX * Core::MAX_SUB);
        float state[Core::STATE_COUNT];
        resetState(state);
        float previousSubstep = state[0];
        int substepCrossings = 0;
        int endpointCrossings = 0;
        for (int observation = 0; observation < 20000; ++observation) {
            const float previousEndpoint = state[0];
            Core::advanceObservation(
                state, schedule.h, schedule.substeps,
                DEFAULT_CURRENT, Core::R_DEFAULT, DEFAULT_ADAPT,
                [&](const float& acceptedMembrane) {
                    if (upwardSpikeCrossing(previousSubstep, acceptedMembrane))
                        ++substepCrossings;
                    previousSubstep = acceptedMembrane;
                });
            if (upwardSpikeCrossing(previousEndpoint, state[0]))
                ++endpointCrossings;
            Core::repair(state);
        }
        const bool substepEventsOk = substepCrossings > 1000
            && endpointCrossings > 0 && substepCrossings > endpointCrossings;
        std::printf("SPIKE accepted-substep crossings=%d, endpoint-only=%d: %s\n",
                    substepCrossings, endpointCrossings,
                    substepEventsOk ? "PASS" : "FAIL");
        if (!substepEventsOk)
            ++failures;
    }

    const float sampleRate = 48000.f;
    int sweepCount = 0;
    float maxRawAbs = 0.f;
    bool rawClampUsed = false;
    for (float current = 0.4f; current <= 4.001f; current += 0.4f) {
        for (float r = Core::R_MIN; r <= Core::R_MAX + 1e-7f; r *= 1.7f) {
            for (float adapt = 1.f; adapt <= 5.001f; adapt += 1.f) {
                for (float oct = -4.f; oct <= 4.001f; oct += 4.f) {
                    float state[Core::STATE_COUNT];
                    resetState(state);
                    const float pitchHz = coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
                    const coalescent::neuron::ScalarSchedule schedule =
                        coalescent::neuron::scalarSchedule<Core>(pitchHz, sampleRate);
                    const long samples = (long) (sampleRate * 0.4f);
                    for (long i = 0; i < samples; ++i) {
                        Core::advanceObservation(
                            state, schedule.h, schedule.substeps, current, r, adapt);
                        for (int j = 0; j < Core::STATE_COUNT; ++j)
                            maxRawAbs = std::max(maxRawAbs, std::fabs(state[j]));
                        const bool rawSafe = rawStateIsSafe(state);
                        const bool withinClamp = rawStateIsWithinClamp(state);
                        rawClampUsed = rawClampUsed || !withinClamp;
                        const bool retained = Core::repair(state);
                        if (!rawSafe || !retained || !rawStateIsWithinClamp(state)) {
                            std::printf("FAIL: escaped I=%.2f r=%.4f s=%.1f oct=%.0f: x=%g y=%g z=%g\n",
                                        current, r, adapt, oct, state[0], state[1], state[2]);
                            ++failures;
                            break;
                        }
                    }
                    ++sweepCount;
                }
            }
        }
    }
    std::printf("STABILITY: %d combinations finite and bounded; max raw |state|=%.2f; raw clamp %s\n",
                sweepCount, maxRawAbs, rawClampUsed ? "was needed" : "was not needed");

    auto measureAudioHz = [&](float oct, float current, float r, float adapt) -> double {
        float state[Core::STATE_COUNT];
        resetState(state);
        const float pitchHz = coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
        const coalescent::neuron::ScalarSchedule schedule =
            coalescent::neuron::scalarSchedule<Core>(pitchHz, sampleRate);
        for (long i = 0; i < (long) (sampleRate * 0.3f); ++i)
            advanceScheduled(state, schedule, current, r, adapt);

        float previous = state[0];
        std::vector<double> intervals;
        double lastCross = -1.0;
        const long maxSamples = (long) (sampleRate * 2.f);
        for (long sample = 0; sample < maxSamples; ++sample) {
            advanceScheduled(state, schedule, current, r, adapt);
            if (upwardSpikeCrossing(previous, state[0])) {
                if (lastCross >= 0.0)
                    intervals.push_back(sample - lastCross);
                lastCross = sample;
            }
            previous = state[0];
        }
        const double periodSamples = medianWithinBurst(intervals);
        return periodSamples > 0.0 ? sampleRate / periodSamples : -1.0;
    };

    for (float oct = -1.f; oct <= 2.001f; oct += 1.f) {
        const double wanted = coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
        const double measured =
            measureAudioHz(oct, DEFAULT_CURRENT, Core::R_DEFAULT, DEFAULT_ADAPT);
        const double cents = measured > 0.0 ? 1200.0 * std::log2(measured / wanted) : 0.0;
        std::printf("PITCH: oct=%+.0f want=%7.2f Hz got=%7.2f Hz (%+6.1f cents)\n",
                    oct, wanted, measured, cents);
        if (measured < 0.0 || std::fabs(cents) > 60.0) {
            std::printf("  FAIL: spike pitch off by more than 60 cents\n");
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
