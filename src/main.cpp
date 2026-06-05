#include <iostream>
#include <random>
#include "absl/base/log_severity.h"
#include "absl/log/globals.h"
#include "ortools/base/init_google.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"
using namespace operations_research::sat;

#include "../include/Graph.h"

#include <string>
#include "google/protobuf/text_format.h"

void print_model(const operations_research::sat::CpModelBuilder& model) {
    std::string model_str;
    // Serialize the underlying proto to a text string
    if (google::protobuf::TextFormat::PrintToString(model.Proto(), &model_str)) {
        std::cout << "--- CP-SAT Model Structure ---" << std::endl;
        std::cout << model_str << std::endl;
    } else {
        std::cerr << "Failed to serialize the model to text format." << std::endl;
    }
}

void shortest_path_example() {
    std::mt19937 gen(35);

    int source_node = 0;
    int target_node = 99;
    auto graph = Graph::randomGraph(gen, 1000, 0.005, 1, 15);
    while (!graph.hasPath(source_node, target_node)) {
        graph.addRandomEdge(gen, 1, 15);
    }

    CpModelBuilder cp_model;
    // Create variables for each node
    operations_research::Domain boolDomain(0, 1);
    std::unordered_map<int, IntVar> nodeVars;
    const operations_research::Domain nodeVarDomain(0, 25);
    std::unordered_map<int, LinearExpr> neighborCounts;
    for (int node: graph.nodes) {
        nodeVars.insert({node, cp_model.NewIntVar(nodeVarDomain).WithName(std::format("node-{}", node))});
        LinearExpr expr = (node == source_node || node == target_node) ? 1 : 0;
        neighborCounts.insert({node, expr});
    }
    // Now for each edge
    std::unordered_map<EdgeInfo, IntVar> edgeVars;
    LinearExpr edgeCost;
    for (std::pair<const EdgeInfo, int> edge: graph.edges) {
        auto ev = cp_model.NewIntVar(boolDomain).WithName(std::format("edge-{},{}({})", edge.first.n1, edge.first.n2, edge.second));
        edgeVars.insert({edge.first, ev});
        edgeCost += ev * edge.second;
        neighborCounts[edge.first.n1] += ev;
        neighborCounts[edge.first.n2] += ev;
    }

    std::unordered_map<int, IntVar> ys;
    for (int node: graph.nodes) {
        ys[node] = cp_model.NewIntVar(boolDomain).WithName(std::format("y-{}", node));
        cp_model.AddEquality(nodeVars[node], neighborCounts[node]);
        cp_model.AddEquality(nodeVars[node], 2*ys[node]);
    }

    cp_model.Minimize(edgeCost);
    const CpSolverResponse response = Solve(cp_model.Build());

    if (response.status() == OPTIMAL) {
        unordered_map<EdgeInfo, int> solution = {};
        for (auto edge: graph.edges) {
            if (SolutionIntegerValue(response, edgeVars[edge.first]) > 0) {
                solution.insert(edge);
                graph.markEdge(edge.first.n1, edge.first.n2);
            }
        }
        //std::cout << graph.toDot() << std::endl;
        std::cout << "Optimal solution found." << std::endl;
        std::cout << "Cost: " << response.objective_value() << std::endl;
    }
}

template <int LATTICE_SIZE, typename RNG>
std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> make_random_error(RNG& gen, double triangle_chance) {
    std::uniform_real_distribution tri_dist(0.0, 1.0);

    std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> result{};
    for (int x = 0; x < LATTICE_SIZE + 1; x++) {
        for (int y = 0; y < LATTICE_SIZE + 1; y++) {
            result[x][y] = 0;
        }
    }

    for (int x = 0; x < LATTICE_SIZE * 2; x++) {
        for (int y = 0; y < LATTICE_SIZE; y++) {
            if (x % 2 == 0) {
                result[x/2][y] ^= 1;
                result[x/2][y+1] ^= 1;
                result[x/2+1][y] ^= 1;
            } else {
                result[x/2+1][y] ^= 1;
                result[x/2][y+1] ^= 1;
                result[x/2+1][y+1] ^= 1;
            }
        }
    }

    return result;
}

template<int LATTICE_SIZE>
std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> solve_lattice(std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error, const int REGION_SIZE) {
    std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> triangle_toggles{};
    for (int y = 0; y < LATTICE_SIZE; y++) {
        for (int x = 0; x < LATTICE_SIZE * 2; x++) {
            triangle_toggles[y][x] = 0;
        }
    }

    int region_increment = REGION_SIZE + 1;

    for (int region_x = 0; region_x < LATTICE_SIZE * 2; region_x += region_increment * 2) {
        for (int region_y = 0; region_y < LATTICE_SIZE; region_y += region_increment) {
            CpModelBuilder cp_model;
            std::array<std::array<BoolVar, LATTICE_SIZE>, LATTICE_SIZE * 2> tri_vars{};
            LinearExpr total_toggles = 0;
            for (int x = region_x; x < LATTICE_SIZE * 2 && x <= region_x + (REGION_SIZE * 2); x++) {
                for (int y = region_y; y < LATTICE_SIZE && y <= region_y + REGION_SIZE; y++) {
                    tri_vars[x][y] = cp_model.NewBoolVar().WithName(std::format("tri-({},{})", x, y));
                    total_toggles += tri_vars[x][y];
                }
            }
            operations_research::Domain errorDomain(0, 7);
            std::array<std::array<IntVar, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> error_counts{};
            std::array<std::array<IntVar, LATTICE_SIZE + 1>, LATTICE_SIZE + 1>  ks{};
            operations_research::Domain modDomain(0, 1);
            std::array<std::array<IntVar, LATTICE_SIZE + 1>, LATTICE_SIZE + 1>  mods{};
            for (int x = region_x / 2; x < LATTICE_SIZE + 1 && x <= region_x + REGION_SIZE; x++) {
                for (int y = region_y; y < LATTICE_SIZE + 1 && y <= region_y + REGION_SIZE; y++) {
                    LinearExpr error = sample_error[y][x]; // Inverted so that inputting the errors from the decoder is 1:1 with how it appears
                    bool is_boundary_point = x == (region_x / 2) || x == (region_x / 2) + REGION_SIZE || y == region_y || y == region_y + REGION_SIZE;

                    std::array<std::pair<int, int>, 6> neighbors = {{
                        {(2 * x) - 1, y - 1},
                        {(2 * x), y - 1},
                        {(2 * x) + 1, y - 1},
                        {(2 * x) - 2, y},
                        {(2 * x) - 1, y},
                        {(2 * x), y},
                    }};

                    for (std::pair<int, int> neighbor: neighbors) {
                        if (neighbor.first >= 0 && neighbor.first < LATTICE_SIZE * 2 && neighbor.second >= 0 && neighbor.second < LATTICE_SIZE) {
                            error += tri_vars[neighbor.first][neighbor.second];
                        }
                    }

                    if (!is_boundary_point) {
                        error_counts[x][y] = cp_model.NewIntVar(errorDomain).WithName(std::format("error-({},{})", x,y));
                        cp_model.AddEquality(error_counts[x][y], error);

                        ks[x][y] = cp_model.NewIntVar(errorDomain).WithName(std::format("k-({},{})", x, y));
                        mods[x][y] = cp_model.NewIntVar(modDomain).WithName(std::format("mod-({},{})", x, y));
                        cp_model.AddEquality(error_counts[x][y], (2 * ks[x][y]) + mods[x][y]);

                        cp_model.AddEquality(mods[x][y], 0);
                    }
                }
            }

            cp_model.Minimize(total_toggles);

            const CpSolverResponse response = Solve(cp_model.Build());

            for (int x = region_x; x < LATTICE_SIZE * 2 && x <= region_x + (REGION_SIZE * 2); x++) {
                for (int y = region_y; y < LATTICE_SIZE && y <= region_y + REGION_SIZE; y++) {
                    triangle_toggles[y][x] = SolutionIntegerValue(response, tri_vars[x][y]);
                }
            }
        }
    }

    return triangle_toggles;
}

int main() {
    absl::InitializeLog();
    absl::SetStderrThreshold(absl::LogSeverity::kInfo);
    std::mt19937 gen(35);

    constexpr int LATTICE_SIZE = 7;
    std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error = {{
        {1, 1, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 1, 0, 0, 0, 1},
        {0, 1, 1, 1, 0, 0, 1, 1},
        {1, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 1, 1},
        {0, 0, 1, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};
    constexpr int REGION_SIZE = 3;

    std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> triangle_toggles = solve_lattice<LATTICE_SIZE>(sample_error, REGION_SIZE);

    for (int y = 0; y < LATTICE_SIZE; y++) {
        for (int x = 0; x < LATTICE_SIZE * 2; x++) {
            std::cout << triangle_toggles[y][x] << " ";
        }
        std::cout << "\n";
        for (int j = 0; j < y+1; j++) {
            std::cout << " ";
        }
    }

    return 0;
}
