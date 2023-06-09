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

#pragma once
#ifndef DELAUNATOR_HPP
#define DELAUNATOR_HPP

#include <limits>
#include <ostream>
#include <vector>

namespace delaunator
{
constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

class Point final
{
  public:
    explicit Point(double x, double y) : m_x(x), m_y(y)
    {
    }

    Point() : m_x(0.0), m_y(0.0)
    {
    }

    const double &x() const
    {
        return m_x;
    }

    const double &y() const
    {
        return m_y;
    }

    double magnitudeSquared() const noexcept
    {
        return m_x * m_x + m_y * m_y;
    }

    static double determinant(const Point &p1, const Point &p2) noexcept
    {
        return p1.m_x * p2.m_y - p1.m_y * p2.m_x;
    }

    static Point vector(const Point &p1, const Point &p2) noexcept
    {
        return Point{p2.m_x - p1.m_x, p2.m_y - p1.m_y};
    }

    static double distanceSquared(const Point &p1, const Point &p2) noexcept
    {
        Point vec = vector(p1, p2);
        return vec.m_x * vec.m_x + vec.m_y * vec.m_y;
    }

    static bool equal(const Point &p1, const Point &p2, double span) noexcept
    {
        const double dist = distanceSquared(p1, p2) / span;

        return (dist < std::numeric_limits<double>::epsilon());
    }

  private:
    double m_x;
    double m_y;
};

inline std::ostream &operator<<(std::ostream &out, const Point &p)
{
    out << p.x() << "/" << p.y();
    return out;
}

class Points final
{
  public:
    using const_iterator = const Point *;

    explicit Points(const std::vector<double> &coords) : m_coords(coords)
    {
    }

    const Point &operator[](const std::size_t offset) const noexcept
    {
        return reinterpret_cast<const Point &>(*(m_coords.data() + (offset * 2)));
    }

    Points::const_iterator begin() const noexcept
    {
        return reinterpret_cast<const Point *>(m_coords.data());
    }

    Points::const_iterator end() const noexcept
    {
        return reinterpret_cast<const Point *>(m_coords.data() + m_coords.size());
    }

    std::size_t size() const noexcept
    {
        return m_coords.size() / 2;
    }

  private:
    const std::vector<double> &m_coords;
};

std::size_t nextHalfEdge(std::size_t e) noexcept;
std::size_t prevHalfEdge(std::size_t e) noexcept;

std::vector<double> getHullCoordinates(const std::vector<std::size_t> &hull_indices, const std::vector<double> &coords)
{
    std::vector<double> hull_coords;
    hull_coords.reserve(2 * hull_indices.size());

    for (const std::size_t hull_index : hull_indices)
    {
        const double x = coords[2 * hull_index];
        const double y = coords[2 * hull_index + 1];
        hull_coords.push_back(x);
        hull_coords.push_back(y);
    }

    return hull_coords;
}

class Delaunator final
{
  public:
    const std::vector<double> &coords;
    Points m_points;

    // 'triangles' stores the indices to the 'X's of the input
    // 'coords'.
    std::vector<std::size_t> triangles;

    // 'half_edges' store indices into 'triangles'.  If half_edges[X] = Y,
    // It says that there's an edge from X to Y where a) X and Y are
    // both indices into triangles and b) X and Y are indices into different
    // triangles in the array.  This allows you to get from a triangle to
    // its adjacent triangle.  If the a triangle edge has no adjacent triangle,
    // its half edge will be INVALID_INDEX.
    std::vector<std::size_t> half_edges;

    std::vector<std::size_t> hull_prev;
    std::vector<std::size_t> hull_next;

    // This contains indexes into the triangles array.
    std::vector<std::size_t> hull_tri;
    std::size_t hull_start;

    Delaunator(const std::vector<double> &in_coords);
    double getHullArea();
    double getTriangleArea();
    std::vector<std::size_t> getHullIndices();
    double edgeLength(const std::size_t e_a) noexcept;
    std::size_t getInteriorPoint(std::size_t e) noexcept;
    std::vector<double> getHullCoordinates();

  private:
    std::vector<std::size_t> m_hash;
    Point m_center;
    std::size_t m_hash_size;
    std::vector<std::size_t> m_edge_stack;

    std::size_t legalize(std::size_t a);
    std::size_t getHashKey(const double x, const double y) const;
    std::size_t addTriangle(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t a, std::size_t b,
                            std::size_t c);
    void link(std::size_t a, std::size_t b);
};

} // namespace delaunator

#endif // DELAUNATOR_HPP