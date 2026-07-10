#pragma once

#include <atomic>
#include <cstdint>

namespace coalescent {

// Single-producer / single-consumer triple buffer for publishing a display frame
// from the audio thread to the UI thread.
//
//   Producer (audio): fill writable(), then publish().
//   Consumer (UI):    consume() once near the start of draw; the returned
//                     reference stays valid (audio won't reclaim that slot) until
//                     the next consume().
//
// Three slots and one atomic word (middle index + a DIRTY bit) give:
//   • lock-free, allocation-free, no UI-thread copy;
//   • a stable UI frame even if audio publishes many times during one slow draw;
//   • latest-frame (not queue) semantics — what a display wants.
//
// Correctness rests on {readIndex, writeIndex, middle&MASK} always being a
// permutation of {0,1,2}: publish swaps write<->middle, consume swaps read<->middle,
// so the audio-write and UI-read slots never alias. Keeping the index and the
// dirty flag in ONE atomic is essential — two independent atomics admit an
// interleaving where the consumer acknowledges one publication but picks up
// another and later swaps back to an older frame.
template <typename T>
class DisplaySnapshot {
    static constexpr uint32_t DIRTY      = 1u << 31;
    static constexpr uint32_t INDEX_MASK = 0x3u;

    // No extra alignment: display payloads are multi-KB (far larger than a cache
    // line, so false sharing between slots is a non-issue), and over-aligning the
    // slot would over-align the whole Module — which Rack allocates with plain new.
    struct Slot { T value{}; };
    Slot slots[3];

    // middle-buffer index (bits 0-1) + DIRTY. Start: UI reads 0, middle 1, audio writes 2.
    std::atomic<uint32_t> middle{1u};
    uint32_t writeIndex = 2u;   // audio-thread-owned
    uint32_t readIndex  = 0u;   // UI-thread-owned

public:
    DisplaySnapshot() = default;
    DisplaySnapshot(const DisplaySnapshot&) = delete;
    DisplaySnapshot& operator=(const DisplaySnapshot&) = delete;

    // Audio thread: fill this completely, then call publish().
    T& writable() noexcept { return slots[writeIndex].value; }

    // Audio thread: hand the filled buffer to the consumer, take the old middle to write next.
    void publish() noexcept {
        uint32_t previous = middle.exchange(writeIndex | DIRTY, std::memory_order_acq_rel);
        writeIndex = previous & INDEX_MASK;
    }

    // UI thread: latest published frame (or the last one consumed, if none new).
    const T& consume() noexcept {
        uint32_t state = middle.load(std::memory_order_acquire);
        if (state & DIRTY) {
            uint32_t previous = middle.exchange(readIndex, std::memory_order_acq_rel);
            readIndex = previous & INDEX_MASK;
        }
        return slots[readIndex].value;
    }
};

} // namespace coalescent
