#pragma once

#include <set>
#include <tuple>
#include <unordered_map>
#include <random>
#include <algorithm>

#include "EdgeInfo.h"
#include "ipm/ipx/basis.h"

using std::set;
using std::tuple;
using std::unordered_map;

struct Graph {
    set<int> nodes;
    unordered_map<EdgeInfo, int> edges;
    set<EdgeInfo> marked_edges;

    void addNode(const int id) {
        nodes.insert(id);
    }

    void addEdge(const int n1, const int n2, const int weight) {
        const EdgeInfo info(n1, n2);
        edges.insert({info, weight});
    }

    void markEdge(const int n1, const int n2) {
        const EdgeInfo info(n1, n2);
        marked_edges.insert(info);
    }

    std::string toDot() const {
        std::ostringstream os;
        os << "graph G {\n  ";
        for (int node : nodes) {
            os << node << "; ";
        }
        os << "\n  ";
        for (auto edge : edges) {
            os << edge.first.n1 << " -- " << edge.first.n2;
            os << " [label=\"" << edge.second << "\"";
            if (marked_edges.contains(edge.first)) {
                os << " color=\"red\"";
            }
            os << "]; ";
        }
        os << "\n}";
        return os.str();
    }

    template <typename RNG>
    void ensureRandomPath(RNG& gen, const int from, const int to, const int weight_lower_bound, const int weight_upper_bound) {
        std::uniform_int_distribution<int> dist(weight_lower_bound, weight_upper_bound);

        set visited = {from};
        int from_temp = from;
        while (from_temp != to) {
            // Pick a random next node that has not been visited
            int next = from_temp;
            while (visited.contains(next)) {
                std::sample(nodes.begin(), nodes.end(), &next, 1, gen);
            }
            if (!edges.contains(EdgeInfo(next, from_temp)) && !edges.contains(EdgeInfo(from_temp, next))) {
                const int weight = dist(gen);
                addEdge(from_temp, next, weight);
            }
            visited.insert(next);
            from_temp = next;
        }
    }

    template <typename RNG>
    static Graph randomGraph(RNG& gen, const int node_count, const double edge_chance, const int weight_lower_bound, const int weight_upper_bound) {
        Graph graph;
        for (int i = 0; i < node_count; ++i) {
            graph.addNode(i);
        }

        std::uniform_real_distribution edge_dist(0.0, 1.0);
        std::uniform_int_distribution weight_dist(weight_lower_bound, weight_upper_bound);
        for (int i = 0; i < node_count; ++i) {
            for (int j = i+1; j < node_count; ++j) {
                if (edge_dist(gen) < edge_chance) {
                    const int weight = weight_dist(gen);
                    graph.addEdge(i, j, weight);
                }
            }
        }

        return graph;
    }

    template <typename RNG>
    void addRandomEdge(RNG& gen, const int weight_lower_bound, const int weight_upper_bound) {
        // Check if the graph is complete
        bool isComplete = true;
        for (const int node1 : nodes) {
            for (const int node2 : nodes) {
                if (node1 == node2) continue;
                if (edges.contains(EdgeInfo(node1, node2))) continue;
                isComplete = false;
            }
        }
        if (isComplete) {return;}

        std::uniform_int_distribution weight_dist(weight_lower_bound, weight_upper_bound);
        int from = 0;
        std::sample(nodes.begin(), nodes.end(), &from, 1, gen);
        int to = 0;
        std::sample(nodes.begin(), nodes.end(), &to, 1, gen);
        while (from == to || edges.contains(EdgeInfo(from, to))) {
            std::sample(nodes.begin(), nodes.end(), &from, 1, gen);
            std::sample(nodes.begin(), nodes.end(), &to, 1, gen);
        }

        int weight = weight_dist(gen);
        addEdge(from, to, weight);
    }

    bool hasPath(int from, int to) {
        return hasPathInternal(from, to, {});
    }

private:
    bool hasPathInternal(int from, int to, set<int> record) {
        if (!nodes.contains(from) || !nodes.contains(to)) return false;
        if (from == to) return true;

        set<int> neighbors {};
        for (int node : nodes) {
            if (edges.contains(EdgeInfo(from, node)) && !record.contains(node)) {
                neighbors.insert(node);
            }
        }

        for (int neighbor : neighbors) {
            set<int> new_record(record);
            new_record.insert(neighbor);
            if (hasPathInternal(neighbor, to, new_record)) {return true;}
        }

        return false;
    }
};