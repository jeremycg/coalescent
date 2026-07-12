#pragma once

#include <algorithm>
#include <cmath>

namespace coalescent {

// Builds a bounded, near-uniform screen-space polyline between two upward
// crossings of `marker`. Sampling and simplification happen incrementally, so a
// crossing only copies an already prepared path instead of doing an unbounded
// retrospective pass on the audio thread.
template <int BINS>
class CompletedPath {
    static_assert(BINS >= 8, "a completed path needs at least eight points");

public:
    struct Point {
        float x = 0.f;
        float y = 0.f;

        Point() = default;
        Point(float x, float y) : x(x), y(y) {}
    };

private:
    Point work[BINS] = {};
    Point latest[BINS] = {};
    bool workImportant[BINS] = {};
    int workCount = 0;
    int latestCount = 0;
    float workLength = 0.f;

    Point beforePrevious;
    Point previous;
    float previousMarker = 0.f;
    float markerAge = 0.f;
    float xScale = 1.f;
    float yScale = 1.f;
    float minMarkerAge = 0.f;
    float maxMarkerAge = 8192.f;
    float minimumSpacing = 0.15f;
    float workSpacing = 0.15f;
    bool havePrevious = false;
    bool haveBeforePrevious = false;
    bool capturing = false;
    bool valid = false;
    bool closed = false;

    float distance(Point a, Point b) const {
        const float dx = std::abs((b.x - a.x) * xScale);
        const float dy = std::abs((b.y - a.y) * yScale);
        const float hi = std::max(dx, dy);
        const float lo = std::min(dx, dy);
        return hi + 0.375f * lo; // fast screen-metric hypot approximation
    }

    void simplify() {
        // Thin protected turns evenly if they ever consume over one third of the
        // buffer, then retain every other ordinary point. This guarantees that
        // one bounded O(BINS) pass releases space even for a long/noisy burst.
        int importantCount = 0;
        for (int i = 1; i < workCount - 1; ++i)
            importantCount += workImportant[i] ? 1 : 0;
        if (importantCount > BINS / 3) {
            const int originalCount = importantCount;
            const int keepCount = BINS / 3;
            int seen = 0, kept = 0;
            for (int i = 1; i < workCount - 1; ++i) {
                if (workImportant[i]) {
                    ++seen;
                    const int keepThroughHere = seen * keepCount / originalCount;
                    if (keepThroughHere > kept)
                        ++kept;
                    else
                        workImportant[i] = false;
                }
            }
        }

        workSpacing *= 2.f;
        int destination = 1;
        bool keepOrdinary = false;
        for (int i = 1; i < workCount - 1; ++i) {
            if (!workImportant[i])
                keepOrdinary = !keepOrdinary;
            if (workImportant[i] || keepOrdinary) {
                work[destination] = work[i];
                workImportant[destination] = workImportant[i];
                ++destination;
            }
        }
        work[destination] = work[workCount - 1];
        workImportant[destination] = workImportant[workCount - 1];
        ++destination;
        workCount = destination;
    }

    void makeRoom() {
        if (workCount >= BINS)
            simplify();
    }

    void appendForced(Point point, bool important) {
        if (workCount > 0 && distance(work[workCount - 1], point) < 1e-9f) {
            workImportant[workCount - 1] = workImportant[workCount - 1] || important;
            return;
        }
        makeRoom();
        if (workCount >= BINS)
            return; // only possible if every stored point was marked important
        if (workCount > 0)
            workLength += distance(work[workCount - 1], point);
        work[workCount] = point;
        workImportant[workCount] = important;
        ++workCount;
    }

    void appendIfMoved(Point point) {
        if (workCount == 0 || distance(work[workCount - 1], point) >= workSpacing)
            appendForced(point, false);
    }

    void appendReversal(Point point) {
        if (workCount == 0) {
            appendForced(point, true);
            return;
        }
        const float separation = distance(work[workCount - 1], point);
        if (separation < 1e-9f || separation >= minimumSpacing * 0.25f)
            appendForced(point, true);
    }

    bool isReversal(Point current) const {
        if (!haveBeforePrevious)
            return false;
        const float ax = previous.x - beforePrevious.x;
        const float ay = previous.y - beforePrevious.y;
        const float bx = current.x - previous.x;
        const float by = current.y - previous.y;
        return ax * bx < -1e-10f || ay * by < -1e-10f;
    }

    bool finish() {
        if (workCount < 3 || !(workLength > 1e-6f))
            return false;

        const float closingGap = distance(work[0], work[workCount - 1]);
        const float averageStep = workLength / (float) (workCount - 1);
        const bool nextClosed = closingGap <= 2.5f * averageStep;
        const int copyCount = nextClosed && closingGap < 0.5f * averageStep
            ? workCount - 1 : workCount;
        std::copy(work, work + copyCount, latest);
        latestCount = copyCount;
        closed = nextClosed;
        valid = true;
        return true;
    }

    void expireLatest() {
        latestCount = 0;
        valid = false;
        closed = false;
    }

    void begin(Point crossing, Point current, float crossingFraction) {
        // Reuse the spatial scale learned by the previous cycle. If the new
        // regime is much shorter, relax it gradually so resolution recovers.
        if (workCount > 0 && workCount < BINS / 2)
            workSpacing = std::max(minimumSpacing, workSpacing * 0.75f);
        workCount = 0;
        workLength = 0.f;
        capturing = true;
        appendForced(crossing, true);
        if (crossingFraction < 1.f)
            appendIfMoved(current);
    }

public:
    void setMetric(float horizontalScale, float verticalScale) {
        xScale = horizontalScale;
        yScale = verticalScale;
    }

    void setMarkerAges(float minimum, float maximum) {
        minMarkerAge = std::max(minimum, 0.f);
        maxMarkerAge = std::max(maximum, minMarkerAge);
    }

    void setMinimumSpacing(float spacing) {
        minimumSpacing = std::max(spacing, 1e-4f);
        workSpacing = workCount == 0 && !capturing
            ? minimumSpacing : std::max(workSpacing, minimumSpacing);
    }

    void reset(bool clearLatest = true) {
        workCount = 0;
        workLength = 0.f;
        previous = {};
        beforePrevious = {};
        previousMarker = 0.f;
        markerAge = 0.f;
        workSpacing = minimumSpacing;
        havePrevious = false;
        haveBeforePrevious = false;
        capturing = false;
        if (clearLatest) {
            std::fill(latest, latest + BINS, Point{});
            latestCount = 0;
            valid = false;
            closed = false;
        }
    }

    // Returns true only when this call promotes a newly completed path.
    // `ageDelta` is dimensionless simulation time, making stale detection
    // independent of sample rate, oversampling and pitch. A discontinuity (hard
    // sync or numerical backstop) aborts only the path under construction; the
    // last valid display path remains available.
    bool push(float x, float y, float marker, float ageDelta,
              bool allowCommit = true, bool discontinuity = false) {
        const Point current{x, y};
        ageDelta = std::max(ageDelta, 0.f);

        if (discontinuity || !havePrevious) {
            workCount = 0;
            workLength = 0.f;
            capturing = false;
            beforePrevious = previous = current;
            previousMarker = marker;
            markerAge = 0.f;
            havePrevious = true;
            haveBeforePrevious = false;
            return false;
        }

        markerAge += ageDelta;
        const bool reversal = isReversal(current);
        const bool upward = previousMarker <= 0.f && marker > 0.f;
        if (upward) {
            const float denominator = marker - previousMarker;
            const float fraction = denominator > 1e-12f
                ? std::max(0.f, std::min(-previousMarker / denominator, 1.f))
                : 1.f;
            const float crossingAge = markerAge - (1.f - fraction) * ageDelta;

            if (crossingAge >= minMarkerAge) {
                const Point crossing{
                    previous.x + (current.x - previous.x) * fraction,
                    previous.y + (current.y - previous.y) * fraction
                };
                bool committed = false;
                if (capturing && crossingAge <= maxMarkerAge) {
                    if (reversal)
                        appendReversal(previous);
                    appendForced(crossing, true);
                    if (allowCommit)
                        committed = finish();
                } else if (crossingAge > maxMarkerAge) {
                    expireLatest();
                }
                begin(crossing, current, fraction);
                markerAge = (1.f - fraction) * ageDelta;
                beforePrevious = previous;
                previous = current;
                previousMarker = marker;
                haveBeforePrevious = true;
                return committed;
            }
        }

        if (markerAge > maxMarkerAge) {
            workCount = 0;
            workLength = 0.f;
            capturing = false;
            expireLatest();
        } else if (capturing) {
            if (reversal)
                appendReversal(previous);
            appendIfMoved(current);
        }
        beforePrevious = previous;
        previous = current;
        previousMarker = marker;
        haveBeforePrevious = true;
        return false;
    }

    const Point* path() const {
        return latest;
    }

    int pathSize() const {
        return latestCount;
    }

    bool hasPath() const {
        return valid;
    }

    bool isClosed() const {
        return closed;
    }
};

// Sample a completed path at constant visible speed. Closed paths wrap; open
// paths make one out-and-back trip per lap so no off-path closing chord is needed.
// Mapping callbacks put x/y into the display metric (normally panel pixels).
template <typename Point, typename MapX, typename MapY>
Point sampleCompletedPath(const Point* points, int count, bool closed, float lap,
                          MapX mapX, MapY mapY) {
    if (count <= 1)
        return points[0];
    const int segments = closed ? count : count - 1;
    float total = 0.f;
    for (int i = 0; i < segments; ++i) {
        const int j = (i + 1) % count;
        total += std::hypot(mapX(points[j].x) - mapX(points[i].x),
                            mapY(points[j].y) - mapY(points[i].y));
    }
    if (!(total > 1e-6f))
        return points[0];

    lap -= std::floor(lap);
    const float direction = closed ? lap
        : (lap < 0.5f ? lap * 2.f : (1.f - lap) * 2.f);
    const float target = direction * total;
    float before = 0.f;
    for (int i = 0; i < segments; ++i) {
        const int j = (i + 1) % count;
        const float segment = std::hypot(mapX(points[j].x) - mapX(points[i].x),
                                         mapY(points[j].y) - mapY(points[i].y));
        if (i == segments - 1 || before + segment >= target) {
            const float f = segment > 1e-6f
                ? std::max(0.f, std::min((target - before) / segment, 1.f)) : 0.f;
            Point result = points[i];
            result.x += (points[j].x - points[i].x) * f;
            result.y += (points[j].y - points[i].y) * f;
            return result;
        }
        before += segment;
    }
    return points[count - 1];
}

} // namespace coalescent
