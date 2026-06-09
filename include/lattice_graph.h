#pragma once

#include <unordered_map>
#include <functional>
#include <queue>

struct int_pair {
    int first;
    int second;

    int_pair(int first, int second) : first(first), second(second) {}

    bool operator==(const int_pair& rhs) const {
        return first == rhs.first && second == rhs.second;
    }

    [[nodiscard]] int_pair clone() const {
        const auto cloned = int_pair(first, second);
        return cloned;
    }

    static bool is_left_closer(int_pair anchor, int_pair left, int_pair right) {
        int left_dist = std::abs(anchor.first - left.first) + std::abs(anchor.second - left.second);
        int right_dist = std::abs(anchor.first - right.first) + std::abs(anchor.second - right.second);
        return left_dist < right_dist;
    }
};

template <>
struct std::hash<int_pair> {
    std::size_t operator()(const int_pair& key) const noexcept {
        const std::size_t h1 = std::hash<int>()(key.first);
        const std::size_t h2 = std::hash<int>()(key.second);

        return h1 ^ (h2 << 1);
    }
};

class lattice_spanning_tree {
public:
    int_pair root;
    std::unordered_map<int_pair, std::tuple<std::shared_ptr<lattice_spanning_tree>, int_pair, int_pair>> children;

    explicit lattice_spanning_tree(int_pair root) : root(root) {
        children = std::unordered_map<int_pair, std::tuple<std::shared_ptr<lattice_spanning_tree>, int_pair, int_pair>>();
    }

    void print_spanning_tree_depth_first(int padding = 0) {
        for (int i = 0; i < padding; i++) {
            std::cout << "  ";
        }
        std::cout << root.first << "," << root.second << std::endl;
        for (tuple<std::shared_ptr<lattice_spanning_tree>, int_pair, int_pair> &entry: children | std::views::values) {
            std::get<0>(entry)->print_spanning_tree_depth_first(padding + 1);
        }
    }
};

class lattice_graph {
    std::unordered_map<int_pair, std::vector<std::tuple<int_pair, int_pair, int_pair>>> adjacencies;
public:
    lattice_graph() = default;

    void add_node(int_pair node) {
        if (!adjacencies.contains(node)) {
            adjacencies[node] = {};
        }
    }

    void add_edge(int_pair node1, int_pair node2, int_pair tri1, int_pair tri2, bool bidirectional = true) {
        add_node(node1);
        add_node(node2);
        adjacencies[node1].push_back({node2, tri1, tri2});
        if (bidirectional) adjacencies[node2].push_back({node1, tri1, tri2});
    }

    void remove_node(int_pair node) {
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

    void remove_edge(int_pair node1, int_pair node2, bool bidirectional = true) {
        for (auto &entry: adjacencies[node1]) {
            if (std::get<0>(entry) == node2) {
                std::erase(adjacencies[node1], entry);
            }
        }
        if (bidirectional) remove_edge(node2, node1, false);
    }

    lattice_graph clone() const {
        lattice_graph cloned;
        for (auto &pair: adjacencies) {
            cloned.add_node(pair.first);
            for (auto &entry: pair.second) {
                cloned.add_edge(pair.first.clone(), std::get<0>(entry).clone(), std::get<1>(entry).clone(), std::get<2>(entry).clone());
            }
        }
        return cloned;
    }

    void print_graph() const {
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

    std::shared_ptr<lattice_spanning_tree> to_spanning_tree(int_pair root) {
        std::queue<int_pair> next;
        next.push(root);
        std::unordered_set<int_pair> already_visited{};
        already_visited.insert(root);
        std::unordered_map<int_pair, std::shared_ptr<lattice_spanning_tree>> constructed_trees{};
        constructed_trees.emplace(root, std::make_shared<lattice_spanning_tree>(root));

        while (!next.empty()) {
            int_pair node = next.front();
            next.pop();
            const std::shared_ptr<lattice_spanning_tree> new_tree = constructed_trees[node];
            for (tuple<int_pair, int_pair, int_pair> &entry: adjacencies[node]) {
                int_pair node2 = std::get<0>(entry);
                int_pair tri1 = std::get<1>(entry);
                int_pair tri2 = std::get<2>(entry);

                if (already_visited.contains(node2)) continue;
                already_visited.insert(node2);
                next.push(node2);

                std::shared_ptr<lattice_spanning_tree> next_tree = std::make_shared<lattice_spanning_tree>(node2);
                new_tree->children.emplace(node2, std::make_tuple(next_tree, tri1, tri2));
                constructed_trees.emplace(node2, next_tree);
            }
        }

        return constructed_trees[root];
    }
};