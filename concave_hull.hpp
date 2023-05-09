// MIT License

// Copyright (c) 2021 Logismos

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
#ifndef CONCAVE_HULL_HPP
#define CONCAVE_HULL_HPP

#include "delaunator.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

namespace geometry
{
template <typename CoordinateType> using ConcaveHullPoint = std::array<CoordinateType, 2>;

template <typename CoordinateType> class ConcaveHull
{
    using Point = ConcaveHullPoint<CoordinateType>;
    using HeapPair = std::pair<std::size_t, double>;

  public:
    explicit ConcaveHull(const std::vector<Point> &points, double chi = 0.1) : chi_(chi)
    {
        if (chi < 0.0 || chi > 1.0)
        {
            throw std::invalid_argument("Chi factor must be between 0 and 1 inclusive");
        }

        coordinates_.reserve(points.size() * 2);
        for (const auto &[x, y] : points)
        {
            coordinates_.push_back(static_cast<double>(x));
            coordinates_.push_back(static_cast<double>(y));
        }

        constructConcaveHull();
    }

    explicit ConcaveHull(const std::vector<CoordinateType> &coords, double chi = 0.1) : chi_(chi)
    {
        if (chi < 0.0 || chi > 1.0)
        {
            throw std::invalid_argument("Chi factor must be between 0 and 1 inclusive");
        }

        coordinates_.reserve(coords.size());
        for (const auto coord : coords)
        {
            coordinates_.push_back(static_cast<double>(coord));
        }

        constructConcaveHull();
    }

    std::vector<std::size_t> getHullIndices() const noexcept
    {
        return hull_indices_;
    }

    std::vector<Point> getHullPoints()
    {
        std::vector<Point> points;
        points.reserve(hull_indices_.size());

        for (const std::size_t hull_index : hull_indices_)
        {
            points.push_back(Point{coordinates_[2 * hull_index], coordinates_[2 * hull_index + 1]});
        }

        return points;
    }

    std::vector<CoordinateType> getHullCoordinates()
    {
        std::vector<CoordinateType> coords;
        coords.reserve(hull_indices_.size() * 2);

        for (const std::size_t hull_index : hull_indices_)
        {
            coords.push_back(coordinates_[2 * hull_index]);
            coords.push_back(coordinates_[2 * hull_index + 1]);
        }

        return coords;
    }

  private:
    static inline bool compare(const HeapPair &left, const HeapPair &right) noexcept
    {
        return (left.second < right.second);
    }

    void constructConcaveHull()
    {
        if (coordinates_.size() < 3)
        {
            return;
        }

        hull_indices_.clear();

        delaunator::Delaunator delaunator(coordinates_);

        // Determine initial points on outside hull
        std::vector<std::size_t> boundary_indices = delaunator.getHullIndices();
        std::set<std::size_t> boundary_set(boundary_indices.begin(), boundary_indices.end());

        // Make max heap of boundary edges with lengths
        std::vector<HeapPair> boundary_heap;
        boundary_heap.reserve(boundary_indices.size());

        auto max_len = std::numeric_limits<double>::lowest();
        auto min_len = std::numeric_limits<double>::max();

        for (const auto boundary_index : boundary_indices)
        {
            const std::size_t e = delaunator.hull_tri[boundary_index];
            const double len = delaunator.edgeLength(e);

            boundary_heap.emplace_back(e, len);
            std::push_heap(boundary_heap.begin(), boundary_heap.end(), compare);

            min_len = std::min(len, min_len);
            max_len = std::max(len, max_len);
        }

        // Determine length parameter
        const double length_param = chi_ * max_len + (1.0 - chi_) * min_len;

        // Iteratively add points to boundary by iterating over the triangles on the hull
        while (!boundary_heap.empty())
        {
            // Get edge with the largest length
            std::pop_heap(boundary_heap.begin(), boundary_heap.end(), compare);
            const auto [e, len] = boundary_heap.back();
            boundary_heap.pop_back();

            // Length of edge too small for our chi factor
            if (len <= length_param)
            {
                break;
            }

            // Find interior point given edge e (a -> b)
            //       e
            //  b <----- a
            //     \   /
            //  e_b \ / e_a
            //       c
            const std::size_t c = delaunator.getInteriorPoint(e);

            // Point already belongs to boundary
            if (boundary_set.count(c))
            {
                continue;
            }

            // Get two edges connected to interior point
            const std::size_t e_n = delaunator::nextHalfEdge(e);

            //  c -> b
            const std::size_t e_b = delaunator.half_edges[e_n];

            //  a -> c
            const std::size_t e_a = delaunator.half_edges[delaunator::nextHalfEdge(e_n)];

            // Add edges to heap
            const double len_a = delaunator.edgeLength(e_a);
            const double len_b = delaunator.edgeLength(e_b);

            boundary_heap.emplace_back(e_a, len_a);
            std::push_heap(boundary_heap.begin(), boundary_heap.end(), compare);

            boundary_heap.emplace_back(e_b, len_b);
            std::push_heap(boundary_heap.begin(), boundary_heap.end(), compare);

            // Update outer hull and connect new edges
            const std::size_t a = delaunator.triangles[e];
            const std::size_t b = delaunator.triangles[e_n];

            delaunator.hull_next[c] = b;
            delaunator.hull_prev[c] = a;
            delaunator.hull_next[a] = c;
            delaunator.hull_prev[b] = c;

            boundary_set.insert(c);
        }

        hull_indices_ = std::move(delaunator.getHullIndices());
    }

    // In
    double chi_;
    std::vector<double> coordinates_;

    // Out
    std::vector<std::size_t> hull_indices_;
};

} // namespace geometry

#endif // CONCAVE_HULL_HPP