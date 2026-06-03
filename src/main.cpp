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
std::array<int, (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1)> make_random_error(RNG& gen, double triangle_chance) {
    std::uniform_real_distribution tri_dist(0.0, 1.0);

    std::array<int, (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1)> result{};
    // First we create a random triangle lattice with i.i.d. triangle placement
    for (int i = 0; i < LATTICE_SIZE * LATTICE_SIZE * 2; i++) {
        if (tri_dist(gen) < triangle_chance) {
            int triangle_x = i % (LATTICE_SIZE * 2);
            int triangle_y = i / (LATTICE_SIZE * 2);

            if (triangle_x % 2 == 0) {
                int top_left_position = (triangle_y * LATTICE_SIZE * 2) + (triangle_x / 2);
                int bottom_position = ((triangle_y + 1) * LATTICE_SIZE * 2) + (triangle_x / 2);
                int top_right_position = (triangle_y * LATTICE_SIZE * 2) + ((triangle_x / 2) + 1);

                result[top_left_position] += 1;
                result[bottom_position] += 1;
                result[top_right_position] += 1;
            } else {
                int top_position = (triangle_y * LATTICE_SIZE * 2) + ((triangle_x / 2) + 1);
                int bottom_left_position = ((triangle_y + 1) * LATTICE_SIZE * 2) + (triangle_x / 2);
                int bottom_right_position = ((triangle_y + 1) * LATTICE_SIZE * 2) + ((triangle_x / 2) + 1);

                result[top_position] += 1;
                result[bottom_left_position] += 1;
                result[bottom_right_position] += 1;
            }
        }
    }

    for (int i = 0; i < (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1); i++) {
        result[i] = result[i] % 2;
    }

    return result;
}

int main() {
    absl::InitializeLog();
    absl::SetStderrThreshold(absl::LogSeverity::kInfo);
    std::mt19937 gen(35);

    constexpr int LATTICE_SIZE = 8;
    std::array<int, (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1)> sample_error = {
        0,0,0,0,0,1,1,0,1,
         0,0,0,1,0,0,1,0,1,
          0,1,0,1,1,1,1,0,1,
           1,1,1,0,1,1,1,1,0,
            1,1,0,0,0,0,1,0,0,
             1,1,1,0,1,1,0,0,0,
              0,0,0,1,1,1,0,1,0,
               1,0,0,0,1,1,1,0,0,
                0,0,1,0,1,1,0,0,0,
    };

    CpModelBuilder cp_model;
    std::array<BoolVar, LATTICE_SIZE * LATTICE_SIZE * 2> tri_vars{};
    LinearExpr total_toggles = 0;
    for (int i = 0; i < LATTICE_SIZE * LATTICE_SIZE * 2; i++) {
        tri_vars[i] = cp_model.NewBoolVar().WithName(std::format("tri-{}", i));
        total_toggles += tri_vars[i];
    }
    operations_research::Domain errorDomain(0, 7);
    std::array<IntVar, (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1)> error_counts{};
    std::array<IntVar, (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1)> ks{};
    operations_research::Domain modDomain(0, 1);
    std::array<IntVar, (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1)> mods{};
    for (int i = 0; i < (LATTICE_SIZE + 1) * (LATTICE_SIZE + 1); i++) {
        error_counts[i] = cp_model.NewIntVar(errorDomain).WithName(std::format("error-{}", i));
        LinearExpr error = sample_error[i];
        int point_x = i % (LATTICE_SIZE + 1);
        int point_y = i / (LATTICE_SIZE + 1);

        int top_left_x = 2 * point_x - 1;
        int top_left_y = point_y - 1;
        if (top_left_x >= 0 && top_left_x < LATTICE_SIZE * 2 && top_left_y >= 0 && top_left_y < LATTICE_SIZE) {
            error += tri_vars[(top_left_y * LATTICE_SIZE * 2) + top_left_x];
        }

        int top_x = 2 * point_x;
        int top_y = point_y - 1;
        if (top_x >= 0 && top_x < LATTICE_SIZE * 2 && top_y >= 0 && top_y < LATTICE_SIZE) {
            error += tri_vars[(top_y * LATTICE_SIZE * 2) + top_x];
        }

        int top_right_x = 2 * point_x + 1;
        int top_right_y = point_y - 1;
        if (top_right_x >= 0 && top_right_x < LATTICE_SIZE * 2 && top_right_y >= 0 && top_right_y < LATTICE_SIZE) {
            error += tri_vars[(top_right_y * LATTICE_SIZE * 2) + top_right_x];
        }

        int bottom_left_x = 2 * point_x - 2;
        int bottom_left_y = point_y;
        if (bottom_left_x >= 0 && bottom_left_x < LATTICE_SIZE * 2 && bottom_left_y >= 0 && bottom_left_y < LATTICE_SIZE) {
            error += tri_vars[(bottom_left_y * LATTICE_SIZE * 2) + bottom_left_x];
        }

        int bottom_x = 2 * point_x - 1;
        int bottom_y = point_y;
        if (bottom_x >= 0 && bottom_x < LATTICE_SIZE * 2 && bottom_y >= 0 && bottom_y < LATTICE_SIZE) {
            error += tri_vars[(bottom_y * LATTICE_SIZE * 2) + bottom_x];
        }

        int bottom_right_x = 2 * point_x;
        int bottom_right_y = point_y;
        if (bottom_right_x >= 0 && bottom_right_x < LATTICE_SIZE * 2 && bottom_right_y >= 0 && bottom_right_y < LATTICE_SIZE) {
            error += tri_vars[(bottom_right_y * LATTICE_SIZE * 2) + bottom_right_x];
        }

        cp_model.AddEquality(error_counts[i], error);

        ks[i] = cp_model.NewIntVar(errorDomain).WithName(std::format("k-{}", i));
        mods[i] = cp_model.NewIntVar(modDomain).WithName(std::format("mod-{}", i));
        cp_model.AddEquality(error_counts[i], (2 * ks[i]) + mods[i]);

        cp_model.AddEquality(mods[i], 0);
    }

    cp_model.Minimize(total_toggles);

    print_model(cp_model);

    const CpSolverResponse response = Solve(cp_model.Build());

    if (response.status() == OPTIMAL) {
        for (int i = 0; i < LATTICE_SIZE * LATTICE_SIZE * 2; i++) {
            std::cout << SolutionIntegerValue(response, tri_vars[i]) << " ";
            if ((i + 1) % (LATTICE_SIZE * 2) == 0) {
                std::cout << "\n";
                for (int j = 0; j < (i + 1) / (LATTICE_SIZE * 2); j++) {
                    std::cout << " ";
                }
            }
        }
        std::cout << std::endl;
        std::cout << "Optimal solution found." << std::endl;
        std::cout << "Cost: " << response.objective_value() << std::endl;
    }

    return 0;
}
