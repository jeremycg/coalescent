#include "../src/dsp/completed_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

static void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "completed_path: %s\n", message);
        std::exit(1);
    }
}

int main() {
    constexpr float PI = 3.14159265358979323846f;

    // A densely and irregularly sampled closed orbit should become a stable,
    // bounded loop without publishing the startup fragment.
    coalescent::CompletedPath<64> circle;
    circle.setMarkerAges(0.01f, 100.f);
    for (int i = -7; i <= 3 * 257 + 9; ++i) {
        const float phase = 2.f * PI * i / 257.f;
        const float warp = phase + 0.015f * std::sin(3.f * phase);
        circle.push(std::cos(warp), std::sin(warp), std::sin(phase), 0.025f);
    }
    require(circle.hasPath(), "a complete second crossing did not publish");
    require(circle.isClosed(), "a closed orbit was classified as open");
    const auto* loop = circle.path();
    const int circleCount = circle.pathSize();
    float minRadius = 10.f, maxRadius = 0.f;
    float maxSegment = 0.f;
    for (int i = 0; i < circleCount; ++i) {
        const float radius = std::sqrt(loop[i].x * loop[i].x + loop[i].y * loop[i].y);
        minRadius = std::min(minRadius, radius);
        maxRadius = std::max(maxRadius, radius);
        const int j = (i + 1) % circleCount;
        const float dx = loop[j].x - loop[i].x;
        const float dy = loop[j].y - loop[i].y;
        const float segment = std::sqrt(dx * dx + dy * dy);
        maxSegment = std::max(maxSegment, segment);
    }
    require(minRadius > 0.998f && maxRadius < 1.002f, "circle sampling left the source polyline");
    require(maxSegment < 0.25f, "closed-loop simplification is too coarse");

    // A one-sample excursion must survive repeated bounded-buffer simplification.
    coalescent::CompletedPath<32> tail;
    tail.setMarkerAges(0.f, 2000.f);
    tail.push(0.f, 0.f, -1.f, 1.f);
    tail.push(0.f, 0.f, 1.f, 1.f);
    for (int i = 1; i <= 1000; ++i) {
        const float x = i == 500 ? 10.f : 0.f;
        tail.push(x, i * 0.01f, i == 1000 ? -1.f : 1.f, 1.f);
    }
    tail.push(0.f, 10.01f, 1.f, 1.f);
    require(tail.hasPath(), "simplified path did not publish");
    float maxTail = -100.f;
    for (int i = 0; i < tail.pathSize(); ++i)
        maxTail = std::max(maxTail, tail.path()[i].x);
    require(maxTail > 9.f, "simplification removed a narrow excursion");

    // Even an adversarial reversal on every sample must stay bounded and leave
    // room for the closing marker (protected turns cannot fill the buffer forever).
    coalescent::CompletedPath<16> zigzag;
    zigzag.setMarkerAges(0.f, 1000.f);
    zigzag.push(0.f, 0.f, -1.f, 1.f);
    zigzag.push(0.f, 0.f, 1.f, 1.f);
    for (int i = 1; i <= 500; ++i)
        zigzag.push(i & 1 ? 1.f : -1.f, i * 0.01f, i == 500 ? -1.f : 1.f, 1.f);
    zigzag.push(0.f, 5.01f, 1.f, 1.f);
    require(zigzag.hasPath() && zigzag.pathSize() <= 16,
            "reversal-heavy path exhausted its bounded buffer");

    // Open paths remain open so the UI can ping-pong instead of drawing a false chord.
    coalescent::CompletedPath<16> open;
    open.setMarkerAges(0.f, 100.f);
    open.push(0.f, 0.f, -1.f, 1.f);
    open.push(0.f, 0.f, 1.f, 1.f);
    for (int i = 1; i <= 10; ++i)
        open.push((float) i, 0.f, i == 10 ? -1.f : 1.f, 1.f);
    open.push(11.f, 0.f, 1.f, 1.f);
    require(open.hasPath() && !open.isClosed(), "an open path was forced closed");

    // A discontinuity keeps the last guide during its grace period but must not
    // connect pre-reset and post-reset states.
    const float oldStart = open.path()[0].x;
    open.push(100.f, 0.f, -1.f, 1.f, true, true);
    require(open.hasPath() && open.path()[0].x == oldStart,
            "discontinuity immediately discarded the completed path");
    open.push(100.f, 0.f, 1.f, 1.f);
    for (int i = 1; i <= 10; ++i)
        open.push(100.f + i, 0.f, i == 10 ? -1.f : 1.f, 1.f);
    open.push(111.f, 0.f, 1.f, 1.f);
    require(open.path()[0].x >= 99.f, "discontinuity introduced a chord from the old state");

    // Suppressed commits discard completed work without invalidating an older path.
    const float heldStart = open.path()[0].x;
    open.push(112.f, 0.f, -1.f, 1.f);
    open.push(113.f, 0.f, 1.f, 1.f, false);
    require(open.path()[0].x == heldStart, "commit throttling changed the published path");

    // With no accepted marker for the configured model-time limit, the guide
    // must expire rather than orbit forever after the live system has reached rest.
    open.push(113.f, 0.f, 1.f, 1000.f);
    require(!open.hasPath() && open.pathSize() == 0, "stale completed path did not expire");

    using Point = coalescent::CompletedPath<16>::Point;
    const Point line[2] = {Point(0.f, 0.f), Point(10.f, 0.f)};
    auto identity = [](float value) { return value; };
    Point head = coalescent::sampleCompletedPath(line, 2, false, 0.25f, identity, identity);
    require(std::fabs(head.x - 5.f) < 1e-5f, "open-path outward sampling is wrong");
    head = coalescent::sampleCompletedPath(line, 2, false, 0.75f, identity, identity);
    require(std::fabs(head.x - 5.f) < 1e-5f, "open-path return sampling is wrong");

    std::puts("completed_path: ok");
    return 0;
}
