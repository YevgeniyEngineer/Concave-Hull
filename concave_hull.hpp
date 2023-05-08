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
std::vector<std::size_t> constructConcaveHull(const std::vector<double> &coords, double chi = 0.1)
{
    if (chi < 0.0 || chi > 1.0)
    {
        throw std::invalid_argument("Chi factor must be between 0 and 1 inclusive");
    }

    delaunator::Delaunator d(coords);

    // Determine initial points on outside hull
    std::vector<std::size_t> boundary_indices = d.getHullIndices();
    std::set<std::size_t> boundary_set(boundary_indices.begin(), boundary_indices.end());

    // Make max heap of boundary edges with lengths
    using HeapPair = std::pair<std::size_t, double>;

    auto cmp = [](HeapPair left, HeapPair right) noexcept { return (left.second < right.second); };

    std::vector<HeapPair> boundary_heap;
    boundary_heap.reserve(boundary_indices.size());

    auto max_len = std::numeric_limits<double>::lowest();
    auto min_len = std::numeric_limits<double>::max();

    for (const auto boundary_index : boundary_indices)
    {
        const std::size_t e = d.hull_tri[boundary_index];
        const double len = d.edgeLength(e);

        boundary_heap.emplace_back(e, len);
        std::push_heap(boundary_heap.begin(), boundary_heap.end(), cmp);

        min_len = std::min(len, min_len);
        max_len = std::max(len, max_len);
    }

    // Determine length parameter
    const double length_param = chi * max_len + (1 - chi) * min_len;

    // Iteratively add points to boundary by iterating over the triangles on the hull
    while (!boundary_heap.empty())
    {
        // Get edge with the largest length
        std::pop_heap(boundary_heap.begin(), boundary_heap.end(), cmp);
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
        const std::size_t c = d.getInteriorPoint(e);

        // Point already belongs to boundary
        if (boundary_set.count(c))
        {
            continue;
        }

        // Get two edges connected to interior point
        const std::size_t e_n = delaunator::nextHalfEdge(e);

        //  c -> b
        const std::size_t e_b = d.half_edges[e_n];

        //  a -> c
        const std::size_t e_a = d.half_edges[delaunator::nextHalfEdge(e_n)];

        // Add edges to heap
        const double len_a = d.edgeLength(e_a);
        const double len_b = d.edgeLength(e_b);

        boundary_heap.emplace_back(e_a, len_a);
        std::push_heap(boundary_heap.begin(), boundary_heap.end(), cmp);

        boundary_heap.emplace_back(e_b, len_b);
        std::push_heap(boundary_heap.begin(), boundary_heap.end(), cmp);

        // Update outer hull and connect new edges
        const std::size_t a = d.triangles[e];
        const std::size_t b = d.triangles[e_n];

        d.hull_next[c] = b;
        d.hull_prev[c] = a;
        d.hull_next[a] = d.hull_prev[b] = c;

        boundary_set.insert(c);
    }

    return d.getHullIndices();
}
} // namespace geometry

#endif // CONCAVE_HULL_HPP