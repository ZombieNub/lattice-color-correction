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

int main() {
    absl::InitializeLog();
    absl::SetStderrThreshold(absl::LogSeverity::kInfo);

    shortest_path_example();

    return 0;
}
