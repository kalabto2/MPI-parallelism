#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

static int max_value = 0, recursion_count = 0;

enum VerticeColor {
    Undefined,
    Red,
    Blue,
};

// Check if it can be valid bipartite graph
bool is_valid(vector<VerticeColor>& vertices,
              vector<pair<int, pair<size_t, size_t>>>& current_edges) {
    for (auto i : current_edges) {
        auto first = vertices.at(i.second.first);
        auto second = vertices.at(i.second.second);
        if (first == Undefined || second == Undefined) continue;
        if (first == second) return false;
    }
    return true;
}

int get_value(vector<pair<int, pair<size_t, size_t>>>& edges) {
    int sum = 0;
    for (auto i : edges) {
        sum += i.first;
    }
    return sum;
}

void DFS(vector<VerticeColor> vertices,
         vector<pair<int, pair<size_t, size_t>>> remainig_edges,
         vector<pair<int, pair<size_t, size_t>>> current_edges);

void print_graph(vector<pair<int, pair<size_t, size_t>>>& edges) {
    for (auto i : edges) {
        cout << "Value: " << i.first << " of (" << i.second.first << ", "
             << i.second.second << ")" << endl;
    }
}

// Run whole program
//  - vertices - saves info about coloring of the edges
//  - remaining_edges - "queue" of sorted edges which consists of <value,
//  <index_of_vertex, index_of_vertex>>
//  - current_edges - "buffer" of selected edges from the remainig_edges
void DFS(vector<VerticeColor> vertices,
         vector<pair<int, pair<size_t, size_t>>> remainig_edges,
         vector<pair<int, pair<size_t, size_t>>> current_edges) {
    // Check if it is valid bipartite graph
    if (!is_valid(vertices, current_edges)) {
        return;
    }

    recursion_count++;

    // We are at the end and found nothing OR it does not make any sense to
    // continue counting
    int current_value = get_value(current_edges),
            remaining_value = get_value(remainig_edges);

    if (remainig_edges.size() == 0 ||
        (current_value + remaining_value) <= max_value) {
        // Ok have we got anything better?
        if (current_value > max_value) {
            max_value = current_value;
        }
        return;
    }

    // Get front and try to color it
    pair<int, pair<size_t, size_t>> front = remainig_edges[0];
    size_t first_vertice = front.second.first,
            second_vertice = front.second.second;

    // remove currently observing edge
    remainig_edges.erase(remainig_edges.begin());

    int counter = 0;
    // First variant
    vector<VerticeColor> new_vertices1(vertices);
    vector<pair<int, pair<size_t, size_t>>> new_remainig_edges1(remainig_edges);
    vector<pair<int, pair<size_t, size_t>>> new_current_edges1(current_edges);
    new_current_edges1.push_back(front);
    new_vertices1.at(first_vertice) = Red;
    new_vertices1.at(second_vertice) = Blue;
    if (is_valid(new_vertices1, new_current_edges1)) {
        DFS(new_vertices1, new_remainig_edges1, new_current_edges1);
        counter++;
    }

    // Second variant
    vector<VerticeColor> new_vertices2(vertices);
    vector<pair<int, pair<size_t, size_t>>> new_remainig_edges2(remainig_edges);
    vector<pair<int, pair<size_t, size_t>>> new_current_edges2(current_edges);
    new_current_edges2.push_back(front);
    new_vertices2.at(first_vertice) = Blue;
    new_vertices2.at(second_vertice) = Red;
    if (is_valid(new_vertices2, new_current_edges2)) {
        DFS(new_vertices2, new_remainig_edges2, new_current_edges2);
        counter++;
    }

    // Third variant
    if (counter != 2) DFS(vertices, remainig_edges, current_edges);
}

int main(void) {
    size_t graph_size, i, z;
    int total_value = 0;
    vector<VerticeColor> vertices;
    // Pairs of value and vertices
    vector<pair<int, pair<size_t, size_t>>> edges;

    // Load input data
    cin >> graph_size;

    vertices.resize(graph_size, VerticeColor::Undefined);
    for (i = 0; i < graph_size; ++i) {
        for (z = 0; z < graph_size; ++z) {
            int value = 0;
            cin >> value;
            if (value > 0 && z > i) {
                edges.emplace_back(make_pair(value, make_pair(i, z)));
                total_value += value;
            }
        }
    }

    cout << "Graph has: " << graph_size << " vertices and: " << edges.size()
         << " edges with total value of: " << total_value
         << "is valid bipartite graph?" << is_valid(vertices, edges) << endl;

    if (edges.size() < 1) {
        assert("invalid input data provided" == nullptr);
        return 1;
    }

    std::sort(std::begin(edges), std::end(edges),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    print_graph(edges);

    vector<pair<int, pair<size_t, size_t>>> tmp;
    DFS(vertices, edges, tmp);

    cout << "Max value of bipartite graph " << max_value
         << " recursion count: " << recursion_count << endl;

    return 0;
}