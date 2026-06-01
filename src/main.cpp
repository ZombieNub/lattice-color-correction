#include <iostream>
#include <random>

#include "../include/Graph.h"

using std::set;
using std::tuple;
using std::unordered_map;

int main() {
    std::mt19937 gen(35);

    auto graph = Graph::randomGraph(gen, 25, 0.1, 1, 10);
    while (!graph.hasPath(0, 24)) {
        graph.addRandomEdge(gen, 1, 10);
    }

    std::string dot = graph.toDot();
    std::cout << dot << std::endl;

    return 0;
}