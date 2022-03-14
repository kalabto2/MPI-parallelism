//
// Created by tomas on 18.02.22.
//

#include "Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <numeric>

int recursion_count = 0;

string readFileIntoString(const string& path) {
    ifstream input_file(path);
    return string((istreambuf_iterator<char>(input_file)), istreambuf_iterator<char>());
}

int retrieve_k(const string & fp) {
    stringstream ss(fp);
    string line;
    while(getline(ss, line, '/')){}
    string file_name = line;

    stringstream iss(file_name);
    while(getline(iss, line, '_')){}
    string k_dirty = line;

    stringstream kss(k_dirty);
    getline(kss, line, '.');

    return stoi(line);
}

Graph::Graph(const string &filepath) {
    this->filepath = filepath;
    string file_content = readFileIntoString(filepath);
    stringstream ss(file_content);
    string line;

    this->e = 0;

    // load lines
    if (!file_content.empty())
    {
        int i = 0;
        while(getline(ss,line,'\n')){
            vector<int> node_edges;
            istringstream iss(line);
            string word;

            // store n
            if (i == 0){
                i ++;
                continue;
            }

            // store edges_weight
            while(iss >> word) {
                int num = stoi (word);
                node_edges.push_back(num);
                if (num != 0)
                    this->e ++;
            }

            this->edges_weight.push_back(node_edges);
            i ++;
        }
    }

    // store k and n
    this->n = (int) this->edges_weight.size();
    this->k = retrieve_k(filepath);
    this->e = this->e / 2;

    for (int i = 0; i < edges_weight.size(); i ++){
        for (int j = 0; j < edges_weight[i].size(); j ++) {
            if (edges_weight[i][j] > 0 && j > i) {
                pair<int, int> edge_vertices;
                edge_vertices.first = i;
                edge_vertices.second = j;
                this->edges_vertices.push_back(edge_vertices);
                this->edges_weight_redux.push_back(edges_weight[i][j]);
            }
        }
    }
}

bool Graph::is_bipartite(vector<BIPART> bipartite_nodes, const vector<vector<bool> >& edges) {
    int node_id = 0;
    for (auto & node_neighbors: edges){
        // check if all neighbors of node are opposite part
        int neighbor_id = 0;
        for (auto neighbor: node_neighbors){
            if (neighbor && (bipartite_nodes[node_id] == bipartite_nodes[neighbor_id]) && (node_id != neighbor_id))
                return false;
            neighbor_id ++;
        }
        node_id ++;
    }

    return true;
}

void Graph::is_coherent_aux(const vector<vector<bool>> &edges, int node_idx, vector<bool> & visited){
    vector<bool> neighbours = edges[node_idx];
    visited[node_idx] = true;

    for (int n_idx = 0; n_idx < neighbours.size(); n_idx ++){
        if (neighbours[n_idx] && !visited[n_idx]){
            is_coherent_aux(edges, n_idx, visited);
        }
    }
}

bool Graph::is_coherent(const vector<vector<bool>> &edges) {
    vector<bool> cluster;
    cluster.push_back(true);
    for (int i = 0; i < edges.size() - 1; i ++)
        cluster.push_back(false);

    is_coherent_aux(edges, 0, cluster);

    // are in cluster all nodes
    return all_of(cluster.begin(), cluster.end(), [](bool v) { return v; });
}

int Graph::get_edges_weight(const vector<vector<bool>> &subgraph_edges) {
    int total_weight = 0;
    for (int i = 0; i < subgraph_edges.size(); i ++)
        for (int j = 0; j < subgraph_edges[i].size(); j ++)
            if (subgraph_edges[i][j])
                total_weight += this->edges_weight[i][j];
    return total_weight / 2;
}

vector<bool> translate (const vector<vector<bool>> & subgraph_edges) {
    vector<bool> result;
    for (int i = 0; i < subgraph_edges.size(); i ++) {
        for (int j = 0; j < subgraph_edges[i].size(); j ++){
            if (j > i) {
                result.push_back(subgraph_edges[i][j]);
            }
        }
    }

    return result;
}

int sum_of_previous (int n) {
    int result = 0;
    for (int i = 0; i <= n; i ++){
        result += i;
    }
    return result;
}

vector<vector<bool>> translate_back  (const vector<bool>& edges) {
    vector<vector<bool>> translated_edges;
    int number_of_nodes = (int) edges.size();
    for (int i = 0; i < number_of_nodes; i ++) {
        vector<bool> neighbors;
        for (int j = 0; j < number_of_nodes; j ++) {
            int idx = i * number_of_nodes - sum_of_previous(i) + j;
            neighbors.push_back(edges[idx]);
        }
        translated_edges.push_back(neighbors);
    }

    return translated_edges;
}

vector<bool> Graph::translate_edges_to_vector () {
    vector<bool> result;
    for (int i = 0; i < this->edges_weight.size(); i ++) {
        for (int j = 0; j < this->edges_weight[i].size(); j ++){
            if (j > i && this->edges_weight[i][j] != 0) {
                result.push_back(true);
            }
        }
    }

    return result;
}

vector<vector<bool>> Graph::translate_vector_to_edges  (const vector<bool>& edges) {
    vector<vector<bool>> translated_edges;

    int idx = 0;
    for (int i = 0; i < this->edges_weight.size(); i ++) {
        vector<bool> neighbors;
        for (int j = 0; j < this->edges_weight[i].size(); j ++){
            if (this->edges_weight[i][j] != 0) {
                if (j > i) {
                    neighbors.push_back(edges[idx]);
                    idx ++;
                }
                else
                    neighbors.push_back(false);
            }
            else
                neighbors.push_back(false);
        }
        translated_edges.push_back(neighbors);
    }

    for (int i = 0; i < translated_edges.size(); i ++) {
        for (int j = 0; j < translated_edges[i].size(); j++) {
            if (j < i)
                translated_edges[i][j] = translated_edges[j][i];
        }
    }

    return translated_edges;
}


bool Graph::has_potential(const vector<bool>& subgraph_edges, int maximum, int weight) {
    int curr_e = subgraph_edges.size();

    int total = 0;
    for (int i = curr_e; i < this->e; i ++) {
        total += this->edges_weight_redux[i];
    }
    return weight + total >= maximum;
}

void Graph::calculate_aux(const vector<bool>& subgraph_edges, const vector<BIPART>& bipartite_nodes, int & maximum,
                          vector<vector<vector<bool>>> &solutions, int node_depth, int weight) {
    recursion_count ++;

    // stop conditions
    if (subgraph_edges.size() >= this->e){
        vector<vector<bool>> translated_edges = this->translate_vector_to_edges(subgraph_edges);
        int current_weight = get_edges_weight(translated_edges);

        // save maximum
        if (current_weight >= maximum &&
        is_coherent(translated_edges) &&
        all_of(bipartite_nodes.begin(), bipartite_nodes.end(), [](BIPART i){return i != NONE;})) {
                if (current_weight > maximum)
                    solutions.clear();
                maximum = current_weight;
//                cout << maximum << endl;
                solutions.push_back(translated_edges);
        }
        return;
    }

    // cut if resting sum isn't higher than current maximum
    if (!has_potential(subgraph_edges, maximum, weight))
        return;

    int i = subgraph_edges.size();
    int f = this->edges_vertices[i].first;
    int s = this->edges_vertices[i].second;

    // 1) take edge
    vector<bool> new_subgraph_edges_with = subgraph_edges;
    new_subgraph_edges_with.push_back(true);
    vector<BIPART> new_nodes = bipartite_nodes;
    if ((bipartite_nodes[f] == ONE && bipartite_nodes[s] == TWO) ||
            (bipartite_nodes[f] == TWO && bipartite_nodes[s] == ONE))
        calculate_aux(new_subgraph_edges_with, new_nodes, maximum, solutions, ++ node_depth, weight + this->edges_weight[f][s]);
    else if (bipartite_nodes[f] == NONE &&
            bipartite_nodes[s] != NONE) {
        new_nodes[f] = (bipartite_nodes[s] == ONE ? TWO : ONE);
        calculate_aux(new_subgraph_edges_with, new_nodes, maximum, solutions, ++ node_depth, weight + this->edges_weight[f][s]);
    }
    else if (bipartite_nodes[f] != NONE &&
            bipartite_nodes[s] == NONE) {
        new_nodes[s] = (bipartite_nodes[f] == ONE ? TWO : ONE);
        calculate_aux(new_subgraph_edges_with, new_nodes, maximum, solutions, ++ node_depth, weight + this->edges_weight[f][s]);
    }
    else if (bipartite_nodes[f] == NONE &&
             bipartite_nodes[s] == NONE) {
        new_nodes[f] = ONE;
        new_nodes[s] = TWO;
        calculate_aux(new_subgraph_edges_with, new_nodes, maximum, solutions, ++ node_depth, weight + this->edges_weight[f][s]);
        new_nodes[s] = ONE;
        new_nodes[f] = TWO;
        calculate_aux(new_subgraph_edges_with, new_nodes, maximum, solutions, ++ node_depth, weight + this->edges_weight[f][s]);
    }

    // 2) dont take edge
    vector<bool> new_subgraph_edges_without = subgraph_edges;
    new_subgraph_edges_without.push_back(false);
    calculate_aux(new_subgraph_edges_without, bipartite_nodes, maximum, solutions, ++ node_depth, weight);

}

pair<int, int> Graph::calculate() {
    int maximum = -1;
    vector<bool> subgraph_edges;
    vector<BIPART> bipartite_nodes;
    bipartite_nodes.assign(this->n, BIPART::NONE);
    bipartite_nodes[0] = BIPART::ONE;
    vector<vector<vector<bool>>> solutions;

    clock_t start;
    start = clock(); // get current time

    // calculate
    calculate_aux(subgraph_edges, bipartite_nodes, maximum, solutions, 0, 0);

    // calculate duration in seconds
    double duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;

//    // print outcome
//    for (int i = 0; i < solutions[0].size(); i ++){
//        for (int j = 0; j < solutions[0][i].size(); j++){
//            cout << solutions[0][i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl << accumulate(this->edges_weight_redux.begin(), this->edges_weight_redux.end(), 0) << endl;
    cout << this->filepath << " -- Solution: " << maximum << " (" << solutions.size() << ") "
    << "(" << duration << " s) (recursion - " << recursion_count << " )" << endl;

//    cout << is_bipartite(, solutions[0]) << endl;
    cout << "is coherent? " << is_coherent(solutions[0]) << endl;

    recursion_count = 0;
}

/*

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
*/