#include <iostream>
#include <random>
#include "absl/base/log_severity.h"
#include "absl/log/globals.h"
#include "ortools/base/init_google.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_checker.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"
using namespace operations_research::sat;

#include "../include/Graph.h"
#include "../include/lattice_graph.h"

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
        std::cout << "Optimal solution found." << std::endl;
        std::cout << "Cost: " << response.objective_value() << std::endl;
    }
}

template <int LATTICE_SIZE, typename RNG>
std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> make_random_error(RNG& gen, double triangle_chance) {
    std::uniform_real_distribution tri_dist(0.0, 1.0);

    std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> result{};
    for (int y = 0; y < LATTICE_SIZE + 1; y++) {
        for (int x = 0; x < LATTICE_SIZE + 1; x++) {
            result[y][x] = 0;
        }
    }

    for (int y = 0; y < LATTICE_SIZE; y++) {
        for (int x = 0; x < LATTICE_SIZE * 2; x++) {
            if (tri_dist(gen) < triangle_chance) {
                if (x % 2 == 0) {
                    result[y][x/2] ^= 1;
                    result[y+1][x/2] ^= 1;
                    result[y][x/2+1] ^= 1;
                } else {
                    result[y][x/2+1] ^= 1;
                    result[y+1][x/2] ^= 1;
                    result[y+1][x/2+1] ^= 1;
                }
            }
        }
    }

    return result;
}

std::array<int_pair, 3> get_adj_points_from_tri(int_pair tri) {
    int x = tri.first;
    int y = tri.second;
    if (x % 2 == 0) {
        return std::array<int_pair, 3>{
            {
                {x/2, y},
                {x/2, y+1},
                {x/2+1, y}
            }};
    } else {
        return std::array<int_pair, 3>{
            {
                {x/2+1, y},
                {x/2, y+1},
                {x/2+1, y+1}
            }};
    }
}

std::array<int_pair, 6> get_adj_tris_from_point(const int_pair point) {
    int x = point.first;
    int y = point.second;
    return std::array<int_pair, 6>{{
        {(2 * x) - 1, y - 1},
        {(2 * x), y - 1},
        {(2 * x) + 1, y - 1},
        {(2 * x) - 2, y},
        {(2 * x) - 1, y},
        {(2 * x), y},
    }};
}

template<int LATTICE_SIZE>
std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> solve_lattice(std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error, const int REGION_SIZE, bool is_border_boundary = true) {
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
            operations_research::Domain errorDomain(0, 7);
            std::array<std::array<IntVar, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> error_counts{};
            std::array<std::array<IntVar, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> ks{};
            operations_research::Domain modDomain(0, 1);
            std::array<std::array<IntVar, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> mods{};
            for (int x = region_x / 2; x < LATTICE_SIZE + 1 && x <= region_x + REGION_SIZE; x++) {
                for (int y = region_y; y < LATTICE_SIZE + 1 && y <= region_y + REGION_SIZE; y++) {
                    LinearExpr error = sample_error[y][x]; // Inverted so that inputting the errors from the decoder is 1:1 with how it appears
                    bool is_boundary_point = false;
                    if (is_border_boundary) {
                        is_boundary_point |= x == (region_x / 2);
                        is_boundary_point |= x == (region_x / 2) + REGION_SIZE;
                        is_boundary_point |= y == region_y;
                        is_boundary_point |= y == region_y + REGION_SIZE;
                    } else {
                        is_boundary_point |= x == (region_x / 2) && x != 0;
                        is_boundary_point |= x == (region_x / 2) + REGION_SIZE && x != LATTICE_SIZE;
                        is_boundary_point |= y == region_y && y != 0;
                        is_boundary_point |= y == region_y + REGION_SIZE && y != LATTICE_SIZE;
                    }

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
                            auto tri_var_instance_possibly_uninitialized = tri_vars[neighbor.first][neighbor.second];
                            if (tri_var_instance_possibly_uninitialized.index() == std::numeric_limits<int32_t>::max() || tri_var_instance_possibly_uninitialized.index() == std::numeric_limits<int32_t>::min()) {
                                tri_vars[neighbor.first][neighbor.second] = cp_model.NewBoolVar().WithName(std::format("tri-({},{})", x, y));
                                total_toggles += tri_vars[neighbor.first][neighbor.second];
                            }
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

            CpModelProto built_model = cp_model.Build();

            const std::string error = ValidateCpModel(built_model);
            if (!error.empty()) {
                std::cerr << "Model is invalid: " << error << std::endl;
            }

            const CpSolverResponse response = Solve(built_model);

            for (int x = region_x; x < LATTICE_SIZE * 2 && x <= region_x + (REGION_SIZE * 2); x++) {
                for (int y = region_y; y < LATTICE_SIZE && y <= region_y + REGION_SIZE; y++) {
                    int answer = SolutionIntegerValue(response, tri_vars[x][y]);
                    triangle_toggles[y][x] = answer;
                }
            }
        }
    }

    return triangle_toggles;
}

template<int LATTICE_SIZE>
bool check_boundary_point(int_pair point, int region_x, int region_y, int REGION_SIZE, bool is_border_boundary = true) {
    if (region_x > 0 && region_y == 0) {
        int dummy = 42;
    }

    int x = point.first;
    int y = point.second;
    bool is_boundary_point = false;
    if (is_border_boundary) {
        is_boundary_point |= x == region_x;
        is_boundary_point |= x == region_x + REGION_SIZE;
        is_boundary_point |= y == region_y;
        is_boundary_point |= y == region_y + REGION_SIZE;
    } else {
        is_boundary_point |= x == region_x && x != 0;
        is_boundary_point |= x == region_x + REGION_SIZE && x != LATTICE_SIZE;
        is_boundary_point |= y == region_y && y != 0;
        is_boundary_point |= y == region_y + REGION_SIZE && y != LATTICE_SIZE;
    }
    return is_boundary_point;
}

template<int LATTICE_SIZE>
std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> solve_lattice_dp(std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error, const int REGION_SIZE, bool is_border_boundary = true) {
    std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> triangle_toggles{};
    for (int y = 0; y < LATTICE_SIZE; y++) {
        for (int x = 0; x < LATTICE_SIZE * 2; x++) {
            triangle_toggles[y][x] = 0;
        }
    }

    int region_increment = REGION_SIZE + 1;

    std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sec{};
    for (int y = 0; y < LATTICE_SIZE + 1; y++) {
        for (int x = 0; x < LATTICE_SIZE + 1; x++) {
            sec[y][x] = sample_error[y][x];
        }
    }

    for (int region_x = 0; region_x < LATTICE_SIZE + 1; region_x += region_increment) {
        for (int region_y = 0; region_y < LATTICE_SIZE + 1; region_y += region_increment) {
            std::unordered_set<int_pair> processed_triangles{};
            for (int x = region_x; x < LATTICE_SIZE + 1 && x < region_x + REGION_SIZE; x++) {
                std::unordered_set<int_pair> tris{};
                std::unordered_map<int_pair, int> parities{};
                for (int y = region_y; y < region_y + REGION_SIZE + 1 && y < LATTICE_SIZE + 1; y++) {
                    if (!check_boundary_point<LATTICE_SIZE>({x, y}, region_x, region_y, REGION_SIZE, is_border_boundary)) {
                        parities[{x, y}] = sec[y][x];
                        auto adj_tris = get_adj_tris_from_point({x, y});
                        for (auto neighbor_tri : adj_tris) {
                            if (processed_triangles.contains(neighbor_tri)) continue;
                            int tx = neighbor_tri.first;
                            int ty = neighbor_tri.second;
                            if (tx >= 0 && tx < LATTICE_SIZE * 2 && ty >= 0 && ty < LATTICE_SIZE) {
                                tris.insert(neighbor_tri);
                            }
                        }
                    }
                }

                bool success = false;
                for (int m = 0; m < tris.size(); m++) {
                    std::vector<bool> bools(tris.size(), false);
                    for (int i = 0; i < m; i++) {
                        bools[i] = true;
                    }
                    do {
                        auto new_parities(parities);
                        std::unordered_set<int_pair> flipped_triangles{};
                        int i = 0;
                        for (auto tri : tris) {
                            if (bools[i]) {
                                flipped_triangles.insert(tri);
                                auto adj_points = get_adj_points_from_tri(tri);
                                for (auto neighbor : adj_points) {
                                    if (new_parities.contains(neighbor)) {
                                        new_parities[neighbor] ^= 1;
                                    }
                                }
                            }
                            i++;
                        }
                        success = true;
                        for (auto v : new_parities | std::views::values) {
                            if (v != 0) {
                                success = false;
                                break;
                            }
                        }
                        if (success) {
                            for (auto tri : flipped_triangles) {
                                triangle_toggles[tri.second][tri.first] ^= 1;
                                for (auto point : get_adj_points_from_tri(tri)) {
                                    sec[point.second][point.first] ^= 1;
                                }
                            }
                            for (auto tri : tris) {
                                processed_triangles.insert(tri);
                            }
                        }
                    } while (std::ranges::prev_permutation(bools).found && success == false);
                    if (success) {
                        break;
                    }
                }
            }
        }
    }

    return triangle_toggles;
}

template<int LATTICE_SIZE>
std::pair<std::vector<int_pair>, int> process_spanning_tree(std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error, const std::shared_ptr<lattice_spanning_tree>& spanning_tree) {
    int current_error = sample_error[spanning_tree->root.second][spanning_tree->root.first];
    std::vector<int_pair> flips{};
    for (std::pair<const int_pair, tuple<std::shared_ptr<lattice_spanning_tree>, int_pair, int_pair> > child: spanning_tree->children) {
        //int_pair target = child.first;
        std::shared_ptr<lattice_spanning_tree> next_tree = std::get<0>(child.second);
        int_pair tri1 = std::get<1>(child.second);
        int_pair tri2 = std::get<2>(child.second);

        auto result = process_spanning_tree<LATTICE_SIZE>(sample_error, next_tree);
        for (auto pair : result.first) {
            flips.push_back(pair);
        }
        if (result.second != 0) {
            flips.push_back(tri1);
            flips.push_back(tri2);
            current_error = current_error ^ result.second;
        }
    }
    return {flips, current_error};
}

template <int LATTICE_SIZE>
std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> apply_triangles(std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error, std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> triangle_choices) {
    std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> result{};
    for (int y = 0; y < LATTICE_SIZE + 1; y++) {
        for (int x = 0; x < LATTICE_SIZE + 1; x++) {
            result[y][x] = sample_error[y][x];
        }
    }

    for (int y = 0; y < LATTICE_SIZE; y++) {
        for (int x = 0; x < LATTICE_SIZE * 2; x++) {
            if (x % 2 == 0) {
                result[y][x/2] ^= triangle_choices[y][x];
                result[y+1][x/2] ^= triangle_choices[y][x];
                result[y][x/2+1] ^= triangle_choices[y][x];
            } else {
                result[y][x/2+1] ^= triangle_choices[y][x];
                result[y+1][x/2] ^= triangle_choices[y][x];
                result[y+1][x/2+1] ^= triangle_choices[y][x];
            }
        }
    }

    return result;
}

template <int LATTICE_SIZE>
std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> solve_residual_error(std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> residual_error, const int REGION_SIZE) {
    std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> triangle_toggles{};
    for (int y = 0; y < LATTICE_SIZE; y++) {
        for (int x = 0; x < LATTICE_SIZE * 2; x++) {
            triangle_toggles[y][x] = 0;
        }
    }

    int_pair centroid{LATTICE_SIZE, LATTICE_SIZE / 2};
    int_pair best_tri = {std::numeric_limits<int>::max()/2, std::numeric_limits<int>::max()/2};

    for (int x = REGION_SIZE * 2; x < LATTICE_SIZE * 2; x += (REGION_SIZE + 1) * 2) {
        for (int y = 0; y < LATTICE_SIZE; y++) {
            if (int_pair::is_left_closer(centroid, {x, y}, best_tri)) {
                best_tri = {x, y};
            }
            if (int_pair::is_left_closer(centroid, {x + 1, y}, best_tri)) {
                best_tri = {x + 1, y};
            }
        }
    }

    for (int y = REGION_SIZE; y < LATTICE_SIZE; y += REGION_SIZE + 1) {
        for (int x = 0; x < LATTICE_SIZE; x++) {
            if (int_pair::is_left_closer(centroid, {x, y}, best_tri)) {
                best_tri = {x, y};
            }
        }
    }

    int_pair best_red{0, 0};
    int_pair best_green{0, 0};
    int_pair best_blue{0, 0};

    int first = best_tri.first;
    int second = best_tri.second;
    if (first % 2 == 0) {
        int_pair top_left = {first/2, second};
        int_pair bottom = {first/2, second+1};
        int_pair top_right = {first/2+1, second};

        int color = (second * 2 + (first / 2)) % 3;
        if (color == 0) {
            best_red = top_left;
            best_green = top_right;
            best_blue = bottom;
        } else if (color == 1) {
            best_red = bottom;
            best_green = top_left;
            best_blue = top_right;
        } else if (color == 2) {
            best_red = top_right;
            best_green = bottom;
            best_blue = top_left;
        }
    } else {
        int_pair top = {first/2 + 1, second};
        int_pair bottom_left = {first/2, second + 1};
        int_pair bottom_right = {first/2 + 1, second + 1};

        int color = (second * 2 + (first / 2)) % 3;
        //std::cout << color << std::endl;
        if (color == 0) {
            best_red = bottom_right;
            best_green = top;
            best_blue = bottom_left;
        } else if (color == 1) {
            best_red = bottom_left;
            best_green = bottom_right;
            best_blue = top;
        } else if (color == 2) {
            best_red = top;
            best_green = bottom_left;
            best_blue = bottom_right;
        }
    }

    //std::cout << best_tri.first << " " << best_tri.second << std::endl;
    //std::cout << best_red.first << " " << best_red.second << std::endl;
    //std::cout << best_green.first << " " << best_green.second << std::endl;
    //std::cout << best_blue.first << " " << best_blue.second << std::endl;

    lattice_graph red_graph;
    lattice_graph green_graph;
    lattice_graph blue_graph;

    // Down-right direction
    for (int x = REGION_SIZE; x < LATTICE_SIZE + 1; x += REGION_SIZE + 1) {
        for (int y = 0; y < LATTICE_SIZE + 1; y++) {
            int opposing_x = x + 1;
            int opposing_y = y + 1;
            if (opposing_x < 0 || opposing_x >= LATTICE_SIZE + 1 || opposing_y < 0 || opposing_y >= LATTICE_SIZE + 1) continue;
            int color = (y * 2 + x) % 3;
            if (color == 0) {
                red_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x), y}, {(2 * x) + 1, y});
            } else if (color == 1) {
                green_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x), y}, {(2 * x) + 1, y});
            } else if (color == 2) {
                blue_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x), y}, {(2 * x) + 1, y});
            }
        }
    }

    // Down direction
    for (int x = REGION_SIZE + 1; x < LATTICE_SIZE + 1; x += REGION_SIZE + 1) {
        for (int y = 0; y < LATTICE_SIZE + 1; y++) {
            int opposing_x = x - 1;
            int opposing_y = y + 2;
            if (opposing_x < 0 || opposing_x >= LATTICE_SIZE + 1 || opposing_y < 0 || opposing_y >= LATTICE_SIZE + 1) continue;
            int color = (y * 2 + x) % 3;
            if (color == 0) {
                red_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x) - 1, y}, {(2 * x) - 2, y + 1});
            } else if (color == 1) {
                green_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x) - 1, y}, {(2 * x) - 2, y + 1});
            } else if (color == 2) {
                blue_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x) - 1, y}, {(2 * x) - 2, y + 1});
            }
        }
    }

    // Down-right direction (Horizontal)
    for (int y = REGION_SIZE; y < LATTICE_SIZE + 1; y += REGION_SIZE + 1) {
        for (int x = 0; x < LATTICE_SIZE + 1; x++) {
            int opposing_x = x + 1;
            int opposing_y = y + 1;
            if (opposing_x < 0 || opposing_x >= LATTICE_SIZE + 1 || opposing_y < 0 || opposing_y >= LATTICE_SIZE + 1) continue;
            int color = (y * 2 + x) % 3;
            if (color == 0) {
                red_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x), y}, {(2 * x) + 1, y});
            } else if (color == 1) {
                green_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x), y}, {(2 * x) + 1, y});
            } else if (color == 2) {
                blue_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x), y}, {(2 * x) + 1, y});
            }
        }
    }

    // Up-right direction
    for (int y = REGION_SIZE + 1; y < LATTICE_SIZE + 1; y += REGION_SIZE + 1) {
        for (int x = 0; x < LATTICE_SIZE + 1; x++) {
            int opposing_x = x + 2;
            int opposing_y = y - 1;
            if (opposing_x < 0 || opposing_x >= LATTICE_SIZE + 1 || opposing_y < 0 || opposing_y >= LATTICE_SIZE + 1) continue;
            int color = (y * 2 + x) % 3;
            if (color == 0) {
                red_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x) + 1, y - 1}, {(2 * x) + 2, y - 1});
            } else if (color == 1) {
                green_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x) + 1, y - 1}, {(2 * x) + 2, y - 1});
            } else if (color == 2) {
                blue_graph.add_edge({x, y}, {opposing_x, opposing_y}, {(2 * x) + 1, y - 1}, {(2 * x) + 2, y - 1});
            }
        }
    }

    std::shared_ptr<lattice_spanning_tree> red_spanning_tree = red_graph.to_spanning_tree(best_red);
    std::shared_ptr<lattice_spanning_tree> green_spanning_tree = green_graph.to_spanning_tree(best_green);
    std::shared_ptr<lattice_spanning_tree> blue_spanning_tree = blue_graph.to_spanning_tree(best_blue);

    std::pair<std::vector<int_pair>, int> red_result = process_spanning_tree<LATTICE_SIZE>(residual_error, red_spanning_tree);
    std::pair<std::vector<int_pair>, int> green_result = process_spanning_tree<LATTICE_SIZE>(residual_error, green_spanning_tree);
    std::pair<std::vector<int_pair>, int> blue_result = process_spanning_tree<LATTICE_SIZE>(residual_error, blue_spanning_tree);

    for (auto flip : red_result.first) {
        triangle_toggles[flip.second][flip.first] ^= 1;
    }
    for (auto flip : green_result.first) {
        triangle_toggles[flip.second][flip.first] ^= 1;
    }
    for (auto flip : blue_result.first) {
        triangle_toggles[flip.second][flip.first] ^= 1;
    }

    if (red_result.second == 1) {
        triangle_toggles[best_tri.second][best_tri.first] ^= 1;
    }

    return triangle_toggles;
}

template <int LATTICE_SIZE>
std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> test_sample(std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error, const int REGION_SIZE) {
    std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> triangle_toggles = solve_lattice_dp<LATTICE_SIZE>(sample_error, REGION_SIZE, true);
    std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> residual_error = apply_triangles<LATTICE_SIZE>(sample_error, triangle_toggles);
    if (REGION_SIZE < LATTICE_SIZE) {
        std::array<std::array<int, LATTICE_SIZE * 2>, LATTICE_SIZE> triangle_toggles_2 = solve_residual_error<LATTICE_SIZE>(residual_error, REGION_SIZE);
        std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> remaining_error = apply_triangles<LATTICE_SIZE>(residual_error, triangle_toggles_2);
        for (int y = 0; y < LATTICE_SIZE; y++) {
            for (int x = 0; x < LATTICE_SIZE * 2; x++) {
                triangle_toggles[y][x] = triangle_toggles[y][x] ^ triangle_toggles_2[y][x];
            }
        }
    }
    return triangle_toggles;
}

template <auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f)
{
    if constexpr (Start < End)
    {
        f(std::integral_constant<decltype(Start), Start>());
        constexpr_for<Start + Inc, End, Inc>(f);
    }
}

int main() {
    absl::InitializeLog();
    absl::SetStderrThreshold(absl::LogSeverity::kInfo);
    std::mt19937 gen(52);
    std::freopen("performance_dp.csv", "w", stdout);

    std::cout << "REGION_SIZE,LATTICE_SIZE,p/100,Elapsed Time (ms)" << "\n";

    template for (constexpr int REGION_SIZE : std::views::iota(2, 6)) {
        template for (constexpr int LATTICE_SIZE : std::views::iota(REGION_SIZE, 101)) {
            for (int p = 0; p <= 100; p += 10) {
                double np = static_cast<double>(p) / static_cast<double>(100);
                auto sample_error = make_random_error<LATTICE_SIZE>(gen, np);
                auto start = std::chrono::steady_clock::now();
                auto result = test_sample<LATTICE_SIZE>(sample_error, REGION_SIZE);
                auto end = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                std::cout << REGION_SIZE << ",";
                std::cout << LATTICE_SIZE << ",";
                std::cout << p << ",";
                std::cout << elapsed.count() << "\n";
            }
        }
    }
    /*
    constexpr int LATTICE_SIZE = 7;
    constexpr int REGION_SIZE = 3;
    std::array<std::array<int, LATTICE_SIZE + 1>, LATTICE_SIZE + 1> sample_error = {
        {
            {0, 0, 0, 0, 0, 0, 1, 1},
            {1, 1, 0, 0, 0, 0, 1, 1},
            {1, 0, 0, 0, 1, 0, 0, 1},
            {0, 1, 0, 1, 1, 0, 1, 0},
            {1, 0, 1, 0, 1, 0, 1, 0},
            {1, 0, 0, 0, 0, 1, 0, 1},
            {1, 0, 0, 0, 0, 0, 1, 1},
            {0, 0, 0, 0, 0, 0, 0, 0},
        }};
    auto result = test_sample<LATTICE_SIZE>(sample_error, REGION_SIZE);

    for (int y = 0; y < LATTICE_SIZE; y++) {
        for (int x = 0; x < LATTICE_SIZE * 2; x++) {
            std::cout << result[y][x] << " ";
        }
        std::cout << std::endl;
        for (int i = 0; i < y + 1; i++) {
            std::cout << " ";
        }
    }
    */
}
