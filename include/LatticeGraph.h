#pragma once

#include <unordered_map>
#include <functional>
#include <queue>

struct IntPair {
    int first;
    int second;

    IntPair(int first, int second) : first(first), second(second) {}

    bool operator==(const IntPair& rhs) const {
        return first == rhs.first && second == rhs.second;
    }

    [[nodiscard]] IntPair clone() const {
        const auto cloned = IntPair(first, second);
        return cloned;
    }
};

template <>
struct std::hash<IntPair> {
    std::size_t operator()(const IntPair& key) const noexcept {
        const std::size_t h1 = std::hash<int>()(key.first);
        const std::size_t h2 = std::hash<int>()(key.second);

        return h1 ^ (h2 << 1);
    }
};

class LatticeSpanningTree {
public:
    IntPair root;
    std::unordered_map<IntPair, std::tuple<std::shared_ptr<LatticeSpanningTree>, IntPair, IntPair>> children;

    explicit LatticeSpanningTree(IntPair root) : root(root) {
        children = std::unordered_map<IntPair, std::tuple<std::shared_ptr<LatticeSpanningTree>, IntPair, IntPair>>();
    }

    void printSpanningTreeDepthFirst(int padding = 0) {
        for (int i = 0; i < padding; i++) {
            std::cout << "  ";
        }
        std::cout << root.first << "," << root.second << std::endl;
        for (tuple<std::shared_ptr<LatticeSpanningTree>, IntPair, IntPair> &entry: children | std::views::values) {
            std::get<0>(entry)->printSpanningTreeDepthFirst(padding + 1);
        }
    }
};

class LatticeGraph {
    std::unordered_map<IntPair, std::vector<std::tuple<IntPair, IntPair, IntPair>>> adjacencies;
public:
    LatticeGraph() = default;

    void addNode(IntPair node) {
        if (!adjacencies.contains(node)) {
            adjacencies[node] = {};
        }
    }

    void addEdge(IntPair node1, IntPair node2, IntPair tri1, IntPair tri2, bool bidirectional = true) {
        addNode(node1);
        addNode(node2);
        adjacencies[node1].push_back({node2, tri1, tri2});
        if (bidirectional) adjacencies[node2].push_back({node1, tri1, tri2});
    }

    void removeNode(IntPair node) {
        if (adjacencies.contains(node)) {
            adjacencies.erase(node);
        }
        for (auto &val: adjacencies | std::views::values) {
            for (auto &entry: val) {
                if (std::get<0>(entry) == node) {
                    std::erase(val, entry);
                }
            }
        }
    }

    void removeEdge(IntPair node1, IntPair node2, bool bidirectional = true) {
        for (auto &entry: adjacencies[node1]) {
            if (std::get<0>(entry) == node2) {
                std::erase(adjacencies[node1], entry);
            }
        }
        if (bidirectional) removeEdge(node2, node1, false);
    }

    LatticeGraph clone() const {
        LatticeGraph cloned;
        for (auto &pair: adjacencies) {
            cloned.addNode(pair.first);
            for (auto &entry: pair.second) {
                cloned.addEdge(pair.first.clone(), std::get<0>(entry).clone(), std::get<1>(entry).clone(), std::get<2>(entry).clone());
            }
        }
        return cloned;
    }

    void printGraph() {
        for (auto &pair: adjacencies) {
            for (auto &entry: pair.second) {
                std::cout << pair.first.first << "," << pair.first.second << " -> ";
                std::cout << std::get<0>(entry).first << "," << std::get<0>(entry).second;
                std::cout << " (" << std::get<1>(entry).first << "," << std::get<1>(entry).second << " | ";
                std::cout << std::get<2>(entry).first << "," << std::get<2>(entry).second << ")";
                std::cout << std::endl;
            }
        }
    }

    std::shared_ptr<LatticeSpanningTree> toSpanningTree(IntPair root) {
        std::queue<IntPair> next;
        next.push(root);
        std::unordered_set<IntPair> alreadyVisited{};
        alreadyVisited.insert(root);
        std::unordered_map<IntPair, std::shared_ptr<LatticeSpanningTree>> constructedTrees{};
        constructedTrees.emplace(root, std::make_shared<LatticeSpanningTree>(root));

        while (!next.empty()) {
            IntPair node = next.front();
            next.pop();
            std::shared_ptr<LatticeSpanningTree> newTree = constructedTrees[node];
            for (tuple<IntPair, IntPair, IntPair> &entry: adjacencies[node]) {
                IntPair node2 = std::get<0>(entry);
                IntPair tri1 = std::get<1>(entry);
                IntPair tri2 = std::get<2>(entry);

                if (alreadyVisited.contains(node2)) continue;
                alreadyVisited.insert(node2);
                next.push(node2);

                std::shared_ptr<LatticeSpanningTree> nextTree = std::make_shared<LatticeSpanningTree>(node2);
                newTree->children.emplace(node2, std::make_tuple(nextTree, tri1, tri2));
                constructedTrees.emplace(node2, nextTree);
            }
        }

        return constructedTrees[root];
    }
};