// MIT License

// Copyright (c) 2018 Volodymyr Bilonenko

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "delaunator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace delaunator
{
//@see
// https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op/33333636#33333636
static inline std::size_t fastModulus(const std::size_t i, const std::size_t c) noexcept
{
    return (i >= c) ? (i % c) : i;
}

// Kahan and Babuska summation, Neumaier variant; accumulates less FP error
static inline double sum(const std::vector<double> &x) noexcept
{
    double sum = x[0];
    double err = 0.0;

    for (std::size_t i = 1; i < x.size(); ++i)
    {
        const double k = x[i];
        const double m = sum + k;

        err += std::fabs(sum) >= std::fabs(k) ? (sum - m + k) : (k - m + sum);
        sum = m;
    }

    return sum + err;
}

static inline double distanceSquared(const double ax, const double ay, const double bx, const double by) noexcept
{
    const double dx = ax - bx;
    const double dy = ay - by;

    return dx * dx + dy * dy;
}

static inline double getCircumRadius(const Point &p1, const Point &p2, const Point &p3) noexcept
{
    const Point d = Point::vector(p1, p2);
    const Point e = Point::vector(p1, p3);

    const double det = Point::determinant(d, e);

    if (std::fabs(det) < std::numeric_limits<double>::epsilon())
    {
        return std::numeric_limits<double>::max();
    }

    const double bl = d.magnitudeSquared();
    const double cl = e.magnitudeSquared();

    const Point radius{(((e.y() * bl - d.y() * cl) * 0.5) / det), (((d.x() * cl - e.x() * bl) * 0.5) / det)};

    return radius.magnitudeSquared();
}

static inline double getCircumRadius(const double ax, const double ay, const double bx, const double by,
                                     const double cx, const double cy) noexcept
{
    const double dx = bx - ax;
    const double dy = by - ay;
    const double ex = cx - ax;
    const double ey = cy - ay;

    const double d = dx * ey - dy * ex;

    // Check if points are collinear
    if (std::fabs(d) <= std::numeric_limits<double>::epsilon())
    {
        return std::numeric_limits<double>::max();
    }

    const double bl = dx * dx + dy * dy;
    const double cl = ex * ex + ey * ey;

    const double x = ((ey * bl - dy * cl) * 0.5) / d;
    const double y = ((dx * cl - ex * bl) * 0.5) / d;

    // Return radius squared
    return ((x * x) + (y * y));
}

static inline bool clockwise(const Point &p0, const Point &p1, const Point &p2) noexcept
{
    // Calculate the vectors from p0 to p1 and from p0 to p2
    const Point v0 = Point::vector(p0, p1);
    const Point v1 = Point::vector(p0, p2);

    // Compute the determinant
    const double det = Point::determinant(v0, v1);

    // Check if the determinant is close to zero (within machine epsilon)
    if (std::fabs(det) <= std::numeric_limits<double>::epsilon())
    {
        return false;
    }

    // Return true if the determinant is negative, which indicates clockwise ordering
    return (det < -std::numeric_limits<double>::epsilon());
}

static inline bool clockwise(const double px, const double py, const double qx, const double qy, const double rx,
                             const double ry) noexcept
{
    const Point p0{px, py};
    const Point p1{qx, qy};
    const Point p2{rx, ry};

    return clockwise(p0, p1, p2);
}

static inline bool counterclockwise(const Point &p0, const Point &p1, const Point &p2) noexcept
{
    // Calculate the vectors from p0 to p1 and from p0 to p2
    const Point v0 = Point::vector(p0, p1);
    const Point v1 = Point::vector(p0, p2);

    // Compute the determinant
    const double det = Point::determinant(v0, v1);

    // Check if the determinant is close to zero (within machine epsilon)
    if (std::fabs(det) <= std::numeric_limits<double>::epsilon())
    {
        return false;
    }

    // Return true if the determinant is positive, which indicates counterclockwise ordering
    return (det > std::numeric_limits<double>::epsilon());
}

static inline bool counterclockwise(const double px, const double py, const double qx, const double qy, const double rx,
                                    const double ry) noexcept
{
    Point p0{px, py};
    Point p1{qx, qy};
    Point p2{rx, ry};

    return counterclockwise(p0, p1, p2);
}

static inline Point getCircumCenter(const double ax, const double ay, const double bx, const double by, const double cx,
                                    const double cy) noexcept
{
    const double dx = bx - ax;
    const double dy = by - ay;
    const double ex = cx - ax;
    const double ey = cy - ay;

    const double d = dx * ey - dy * ex;

    if (std::fabs(d) <= std::numeric_limits<double>::epsilon())
    {
        return Point{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    }

    const double bl = dx * dx + dy * dy;
    const double cl = ex * ex + ey * ey;

    const double x = ax + ((ey * bl - dy * cl) * 0.5) / d;
    const double y = ay + ((dx * cl - ex * bl) * 0.5) / d;

    return Point{x, y};
}

static inline bool isInsideCircumCircle(const double ax, const double ay, const double bx, const double by,
                                        const double cx, const double cy, const double px, const double py) noexcept
{
    const double dx = ax - px;
    const double dy = ay - py;
    const double ex = bx - px;
    const double ey = by - py;
    const double fx = cx - px;
    const double fy = cy - py;

    const double ap = dx * dx + dy * dy;
    const double bp = ex * ex + ey * ey;
    const double cp = fx * fx + fy * fy;

    return ((dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx)) <
            -std::numeric_limits<double>::epsilon());
}

constexpr static double EPSILON = std::numeric_limits<double>::epsilon();

static inline bool checkPointsEqual(const double x1, const double y1, const double x2, const double y2) noexcept
{
    return (std::fabs(x1 - x2) <= EPSILON) && (std::fabs(y1 - y2) <= EPSILON);
}

// monotonically increases with real angle, but doesn't need expensive trigonometry
static inline double pseudoAngle(const double dx, const double dy) noexcept
{
    const double p = dx / (std::fabs(dx) + std::fabs(dy));
    return (((dy > 0.0) ? (3.0 - p) : (1.0 + p)) / 4.0); // [0..1)
}

Delaunator::Delaunator(const std::vector<double> &in_coords) : coords(in_coords), m_points(in_coords)
{
    std::size_t n = coords.size() >> 1;

    std::vector<std::size_t> ids(n);
    std::iota(ids.begin(), ids.end(), 0);

    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();

    for (const Point &p : m_points)
    {
        min_x = std::min(p.x(), min_x);
        min_y = std::min(p.y(), min_y);
        max_x = std::max(p.x(), max_x);
        max_y = std::max(p.y(), max_y);
    }

    const double width = max_x - min_x;
    const double height = max_y - min_y;
    const double span = width * width + height * height; // Everything is square dist.

    Point center((min_x + max_x) / 2.0, (min_y + max_y) / 2.0);

    std::size_t i0 = INVALID_INDEX;
    std::size_t i1 = INVALID_INDEX;
    std::size_t i2 = INVALID_INDEX;

    // pick a seed point close to the centroid
    double min_dist = std::numeric_limits<double>::max();

    for (std::size_t i = 0; i < m_points.size(); ++i)
    {
        const Point &p = m_points[i];

        const double d = Point::distanceSquared(center, p);

        if (d < min_dist)
        {
            i0 = i;
            min_dist = d;
        }
    }

    const Point &p0 = m_points[i0];

    min_dist = std::numeric_limits<double>::max();

    // find the point closest to the seed
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i == i0)
        {
            continue;
        }

        const double d = Point::distanceSquared(p0, m_points[i]);

        if ((d < min_dist) && (d > 0.0))
        {
            i1 = i;
            min_dist = d;
        }
    }

    const Point &p1 = m_points[i1];

    double min_radius = std::numeric_limits<double>::max();

    // find the third point which forms the smallest circumcircle
    // with the first two
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i == i0 || i == i1)
        {
            continue;
        }

        const double r = getCircumRadius(p0, p1, m_points[i]);

        if (r < min_radius)
        {
            i2 = i;
            min_radius = r;
        }
    }

    if (!(min_radius < std::numeric_limits<double>::max()))
    {
        throw std::runtime_error("not triangulation");
    }

    const Point &p2 = m_points[i2];

    if (counterclockwise(p0, p1, p2))
    {
        std::swap(i1, i2);
    }

    const double i0x = p0.x();
    const double i0y = p0.y();
    const double i1x = m_points[i1].x();
    const double i1y = m_points[i1].y();
    const double i2x = m_points[i2].x();
    const double i2y = m_points[i2].y();

    m_center = getCircumCenter(i0x, i0y, i1x, i1y, i2x, i2y);

    // Calculate the distances from the center once to avoid having to
    // calculate for each compare.  This used to be done in the comparator,
    // but GCC 7.5+ would copy the comparator to iterators used in the
    // sort, and this was excruciatingly slow when there were many points
    // because you had to copy the vector of distances.
    std::vector<double> dists;
    dists.reserve(m_points.size());
    for (const Point &p : m_points)
    {
        dists.push_back(distanceSquared(p.x(), p.y(), m_center.x(), m_center.y()));
    }

    // sort the points by distance from the seed triangle getCircumCenter
    std::sort(ids.begin(), ids.end(),
              [&dists](const std::size_t i, const std::size_t j) noexcept { return (dists[i] < dists[j]); });

    // initialize a hash table for storing edges of the advancing convex hull
    m_hash_size = static_cast<std::size_t>(std::ceil(std::sqrt(n)));
    m_hash.resize(m_hash_size);
    std::fill(m_hash.begin(), m_hash.end(), INVALID_INDEX);

    // initialize arrays for tracking the edges of the advancing convex hull
    hull_prev.resize(n);
    hull_next.resize(n);
    hull_tri.resize(n);

    hull_start = i0;

    std::size_t hull_size = 3;

    hull_next[i0] = hull_prev[i2] = i1;
    hull_next[i1] = hull_prev[i0] = i2;
    hull_next[i2] = hull_prev[i1] = i0;

    hull_tri[i0] = 0;
    hull_tri[i1] = 1;
    hull_tri[i2] = 2;

    m_hash[getHashKey(i0x, i0y)] = i0;
    m_hash[getHashKey(i1x, i1y)] = i1;
    m_hash[getHashKey(i2x, i2y)] = i2;

    const std::size_t max_triangles = (n < 3) ? 1 : (2 * n - 5);
    triangles.reserve(max_triangles * 3);
    half_edges.reserve(max_triangles * 3);

    addTriangle(i0, i1, i2, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);

    double xp = std::numeric_limits<double>::quiet_NaN();
    double yp = std::numeric_limits<double>::quiet_NaN();

    // Go through points based on distance from the center.
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t i = ids[k];
        const double x = coords[2 * i];
        const double y = coords[2 * i + 1];

        // skip near-duplicate points
        if (k > 0 && checkPointsEqual(x, y, xp, yp))
        {
            continue;
        }

        xp = x;
        yp = y;

        // skip seed triangle points
        if (checkPointsEqual(x, y, i0x, i0y) || checkPointsEqual(x, y, i1x, i1y) || checkPointsEqual(x, y, i2x, i2y))
        {
            continue;
        }

        // find a visible edge on the convex hull using edge hash
        std::size_t start = 0;

        std::size_t key = getHashKey(x, y);
        for (std::size_t j = 0; j < m_hash_size; ++j)
        {
            start = m_hash[fastModulus(key + j, m_hash_size)];
            if (start != INVALID_INDEX && start != hull_next[start])
            {
                break;
            }
        }

        start = hull_prev[start];
        std::size_t e = start;
        std::size_t q;

        // Advance until we find a place in the hull where our current point
        // can be added.
        while (true)
        {
            q = hull_next[e];

            if (Point::equal(m_points[i], m_points[e], span) || Point::equal(m_points[i], m_points[q], span))
            {
                e = INVALID_INDEX;
                break;
            }

            if (counterclockwise(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1]))
            {
                break;
            }

            e = q;

            if (e == start)
            {
                e = INVALID_INDEX;
                break;
            }
        }

        if (e == INVALID_INDEX) // likely a near-duplicate point; skip it
        {
            continue;
        }

        // add the first triangle from the point
        std::size_t t = addTriangle(e, i, hull_next[e], INVALID_INDEX, INVALID_INDEX, hull_tri[e]);

        hull_tri[i] = legalize(t + 2); // Legalize the triangle we just added.
        hull_tri[e] = t;
        ++hull_size;

        // walk forward through the hull, adding more triangles and
        // flipping recursively
        std::size_t next = hull_next[e];
        while (true)
        {
            q = hull_next[next];
            if (!counterclockwise(x, y, coords[2 * next], coords[2 * next + 1], coords[2 * q], coords[2 * q + 1]))
            {
                break;
            }
            t = addTriangle(next, i, q, hull_tri[i], INVALID_INDEX, hull_tri[next]);
            hull_tri[i] = legalize(t + 2);
            hull_next[next] = next; // mark as removed
            --hull_size;
            next = q;
        }

        // walk backward from the other side, adding more triangles and flipping
        if (e == start)
        {
            while (true)
            {
                q = hull_prev[e];
                if (!counterclockwise(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1]))
                {
                    break;
                }
                t = addTriangle(q, i, e, INVALID_INDEX, hull_tri[e], hull_tri[q]);
                legalize(t + 2);
                hull_tri[q] = t;
                hull_next[e] = e; // mark as removed
                --hull_size;
                e = q;
            }
        }

        // update the hull indices
        hull_prev[i] = e;
        hull_start = e;
        hull_prev[next] = i;
        hull_next[e] = i;
        hull_next[i] = next;

        m_hash[getHashKey(x, y)] = i;
        m_hash[getHashKey(coords[2 * e], coords[2 * e + 1])] = e;
    }
}

double Delaunator::getHullArea()
{
    std::vector<double> hull_area;
    std::size_t e = hull_start;
    std::size_t cnt = 1;
    do
    {
        hull_area.push_back((coords[2 * e] - coords[2 * hull_prev[e]]) *
                            (coords[2 * e + 1] + coords[2 * hull_prev[e] + 1]));
        ++cnt;
        e = hull_next[e];
    } while (e != hull_start);
    return sum(hull_area);
}

double Delaunator::getTriangleArea()
{
    std::vector<double> vals;
    vals.reserve(triangles.size() / 3);
    for (std::size_t i = 0; i < triangles.size(); i += 3)
    {
        const double ax = coords[2 * triangles[i]];
        const double ay = coords[2 * triangles[i] + 1];
        const double bx = coords[2 * triangles[i + 1]];
        const double by = coords[2 * triangles[i + 1] + 1];
        const double cx = coords[2 * triangles[i + 2]];
        const double cy = coords[2 * triangles[i + 2] + 1];
        const double val = std::fabs((by - ay) * (cx - bx) - (bx - ax) * (cy - by));
        vals.push_back(val);
    }
    return sum(vals);
}

std::size_t Delaunator::legalize(std::size_t a)
{
    std::size_t i = 0;
    std::size_t ar = 0;
    m_edge_stack.clear();

    // recursion eliminated with a fixed-size stack
    while (true)
    {
        const std::size_t b = half_edges[a];

        /* if the pair of triangles doesn't satisfy the Delaunay condition
         * (p1 is inside the circum circle of [p0, pl, pr]), flip them,
         * then do the same check/flip recursively for the new pair of triangles
         *
         *           pl                    pl
         *          /||\                  /  \
         *       al/ || \bl            al/    \a
         *        /  ||  \              /      \
         *       /  a||b  \    flip    /___ar___\
         *     p0\   ||   /p1   =>   p0\---bl---/p1
         *        \  ||  /              \      /
         *       ar\ || /br             b\    /br
         *          \||/                  \  /
         *           pr                    pr
         */
        const std::size_t a0 = 3 * (a / 3);
        ar = a0 + (a + 2) % 3;

        if (b == INVALID_INDEX)
        {
            if (i > 0)
            {
                --i;
                a = m_edge_stack[i];
                continue;
            }
            else
            {
                // i = INVALID_INDEX;
                break;
            }
        }

        const std::size_t b0 = 3 * (b / 3);
        const std::size_t al = a0 + (a + 1) % 3;
        const std::size_t bl = b0 + (b + 2) % 3;

        const std::size_t p0 = triangles[ar];
        const std::size_t pr = triangles[a];
        const std::size_t pl = triangles[al];
        const std::size_t p1 = triangles[bl];

        const bool illegal =
            isInsideCircumCircle(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr], coords[2 * pr + 1], coords[2 * pl],
                                 coords[2 * pl + 1], coords[2 * p1], coords[2 * p1 + 1]);

        if (illegal)
        {
            triangles[a] = p1;
            triangles[b] = p0;

            auto hbl = half_edges[bl];

            // Edge swapped on the other side of the hull (rare).
            // Fix the halfedge reference
            if (hbl == INVALID_INDEX)
            {
                std::size_t e = hull_start;
                do
                {
                    if (hull_tri[e] == bl)
                    {
                        hull_tri[e] = a;
                        break;
                    }
                    e = hull_prev[e];
                } while (e != hull_start);
            }

            link(a, hbl);
            link(b, half_edges[ar]);
            link(ar, bl);

            const std::size_t br = b0 + (b + 1) % 3;

            if (i < m_edge_stack.size())
            {
                m_edge_stack[i] = br;
            }
            else
            {
                m_edge_stack.push_back(br);
            }
            ++i;
        }
        else
        {
            if (i > 0)
            {
                --i;
                a = m_edge_stack[i];
                continue;
            }
            else
            {
                break;
            }
        }
    }
    return ar;
}

std::size_t Delaunator::getHashKey(const double x, const double y) const
{
    const double dx = x - m_center.x();
    const double dy = y - m_center.y();
    return fastModulus(
        static_cast<std::size_t>(std::llround(std::floor(pseudoAngle(dx, dy) * static_cast<double>(m_hash_size)))),
        m_hash_size);
}

std::size_t Delaunator::addTriangle(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t a, std::size_t b,
                                    std::size_t c)
{
    std::size_t t = triangles.size();
    triangles.push_back(i0);
    triangles.push_back(i1);
    triangles.push_back(i2);
    link(t, a);
    link(t + 1, b);
    link(t + 2, c);
    return t;
}

void Delaunator::link(const std::size_t a, const std::size_t b)
{
    std::size_t s = half_edges.size();
    if (a == s)
    {
        half_edges.push_back(b);
    }
    else if (a < s)
    {
        half_edges[a] = b;
    }
    else
    {
        throw std::runtime_error("Cannot link edge");
    }
    if (b != INVALID_INDEX)
    {
        std::size_t s2 = half_edges.size();
        if (b == s2)
        {
            half_edges.push_back(a);
        }
        else if (b < s2)
        {
            half_edges[b] = a;
        }
        else
        {
            throw std::runtime_error("Cannot link edge");
        }
    }
}

std::size_t nextHalfEdge(std::size_t e) noexcept
{
    return (e % 3 == 2) ? (e - 2) : (e + 1);
}

std::size_t prevHalfEdge(std::size_t e) noexcept
{
    return (e % 3 == 0) ? (e + 2) : (e - 1);
}

std::vector<std::size_t> Delaunator::getHullIndices()
{
    std::vector<std::size_t> hull_pts;

    std::size_t point = hull_start;
    do
    {
        hull_pts.push_back(point);
        point = hull_next[point];
    } while (point != hull_start);

    hull_pts.push_back(hull_start);

    return hull_pts;
}

double Delaunator::edgeLength(const std::size_t e_a) noexcept
{
    const std::size_t e_b = nextHalfEdge(e_a);

    const double x_a = coords[2 * triangles[e_a]];
    const double y_a = coords[2 * triangles[e_a] + 1];

    const double x_b = coords[2 * triangles[e_b]];
    const double y_b = coords[2 * triangles[e_b] + 1];

    return std::sqrt(std::pow(x_a - x_b, 2.0) + std::pow(y_a - y_b, 2.0));
}

std::size_t Delaunator::getInteriorPoint(std::size_t e) noexcept
{
    return triangles[nextHalfEdge(nextHalfEdge(e))];
}

std::vector<double> Delaunator::getHullCoordinates()
{
    const std::vector<std::size_t> hull_indices = getHullIndices();

    std::vector<double> hull_coords;
    hull_coords.reserve(2 * hull_indices.size());

    for (std::size_t hull_index : hull_indices)
    {
        const double x = coords[2 * hull_index];
        const double y = coords[2 * hull_index + 1];
        hull_coords.push_back(x);
        hull_coords.push_back(y);
    }

    return hull_coords;
}

} // namespace delaunator