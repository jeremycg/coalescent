#pragma once

#include "playback.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace coalescent {
namespace lineages {

static const int STATE_VERSION = 1;
static const int FACTORY_LEAF_COUNT = 8;
static const double FACTORY_MUTATE = 0.4;
static const std::uint64_t FACTORY_SEED = UINT64_C(42);
static const std::uint64_t FACTORY_STREAM = UINT64_C(54);

class Hex64Codec {
public:
    enum {
        DIGITS = 16,
        TEXT_SIZE = DIGITS + 1
    };

    static void format(std::uint64_t value, char (&text)[TEXT_SIZE]) {
        static const char digits[] = "0123456789abcdef";
        for (int index = DIGITS - 1; index >= 0; --index) {
            text[index] = digits[value & UINT64_C(0xf)];
            value >>= 4u;
        }
        text[DIGITS] = '\0';
    }

    // Parsing is transactional: destination is untouched unless all sixteen
    // digits are present and valid. Both cases are accepted for compatibility
    // with existing Coalescent patch state; formatting is canonical lowercase.
    static bool parse(const char* text, std::size_t length,
                      std::uint64_t& destination) {
        if (!text || length != static_cast<std::size_t>(DIGITS))
            return false;

        std::uint64_t parsed = 0u;
        for (int index = 0; index < DIGITS; ++index) {
            const char character = text[index];
            unsigned digit = 0u;
            if (character >= '0' && character <= '9')
                digit = static_cast<unsigned>(character - '0');
            else if (character >= 'a' && character <= 'f')
                digit = static_cast<unsigned>(character - 'a' + 10);
            else if (character >= 'A' && character <= 'F')
                digit = static_cast<unsigned>(character - 'A' + 10);
            else
                return false;
            parsed = (parsed << 4u) | static_cast<std::uint64_t>(digit);
        }

        destination = parsed;
        return true;
    }

    static bool parseCString(const char* text, std::uint64_t& destination) {
        if (!text)
            return false;
        std::size_t length = 0u;
        while (length <= static_cast<std::size_t>(DIGITS) && text[length] != '\0')
            ++length;
        return parse(text, length, destination);
    }
};

// Canonical save payload shared by the audio-thread snapshot and JSON restore.
// The Tree carries the generated structural settings. Derived tree fields and
// dirty status are intentionally omitted from the serialized representation:
// validateAndDerive() rebuilds the former and the Rack wrapper recomputes the
// latter from its current SAMPLES/MUTATE controls.
struct Snapshot {
    int version;
    GeneratedState generated;
    double cursor;
    Direction direction;
    bool loop;
    bool running;
    float nodePulseRemaining;
    float mutationPulseRemaining;
    float mrcaPulseRemaining;

    Snapshot()
        : version(STATE_VERSION), cursor(1.0), direction(Direction::Descent),
          loop(true), running(true), nodePulseRemaining(0.f),
          mutationPulseRemaining(0.f), mrcaPulseRemaining(0.f) {}

    static float maximumPulseRemaining() { return 1e-3f; }

    bool validateAndDerive() {
        Snapshot candidate = *this;
        if (candidate.version != STATE_VERSION ||
            !candidate.generated.validateAndDerive() ||
            !std::isfinite(candidate.cursor) ||
            candidate.cursor < 0.0 || candidate.cursor > 1.0 ||
            !validDirection(candidate.direction) ||
            !validPulse(candidate.nodePulseRemaining) ||
            !validPulse(candidate.mutationPulseRemaining) ||
            !validPulse(candidate.mrcaPulseRemaining))
            return false;
        *this = candidate;
        return true;
    }

    // Installs the composition transactionally. Pulse remainders remain in
    // this Snapshot for the Rack wrapper to apply after this method succeeds.
    bool install(KingmanGenerator& generator, Tree& tree,
                 Playback& playback) const {
        Snapshot candidate = *this;
        if (!candidate.validateAndDerive())
            return false;

        KingmanGenerator candidateGenerator = generator;
        Tree candidateTree;
        if (!candidateGenerator.restore(candidate.generated, candidateTree))
            return false;
        Playback candidatePlayback = playback;
        candidatePlayback.restore(candidateTree, candidate.cursor,
                                  candidate.direction, candidate.loop,
                                  candidate.running);

        generator = candidateGenerator;
        tree = candidateTree;
        playback = candidatePlayback;
        return true;
    }

private:
    static bool validDirection(Direction value) {
        return value == Direction::Ancestry || value == Direction::Descent;
    }

    static bool validPulse(float remaining) {
        return std::isfinite(remaining) && remaining >= 0.f &&
               remaining <= maximumPulseRemaining();
    }
};

inline bool makeFactorySnapshot(Snapshot& destination) {
    KingmanGenerator generator(FACTORY_SEED, FACTORY_STREAM);
    Tree tree;
    if (!generator.generate(FACTORY_LEAF_COUNT, FACTORY_MUTATE, tree))
        return false;

    Playback playback;
    playback.setDirection(Direction::Descent);
    playback.setLoop(true);
    playback.installTreeAtSource(tree);

    Snapshot candidate;
    candidate.generated = generator.capture(tree);
    candidate.cursor = playback.cursor();
    candidate.direction = playback.direction();
    candidate.loop = playback.loop();
    candidate.running = playback.running();
    if (!candidate.validateAndDerive())
        return false;
    destination = candidate;
    return true;
}

} // namespace lineages
} // namespace coalescent
