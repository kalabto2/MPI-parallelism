//
// Created by tomas on 18.02.22.
//

#include "Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>

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

//void Graph::is_bipartite_aux(vector<BIPART>& nodes, int node_id) {
//
//}
//
//bool Graph::is_bipartite(const vector<vector<bool> >& edges, int node_number) {
//    vector<BIPART> nodes;
//    nodes.assign(node_number, NONE);
//
//    nodes[0] = ONE;
//    for (int j = 0; j < node_number; j ++) {
//        if (edges[0][j]){
//            if (nodes[0] == ONE)
//                nodes[j] = TWO;
//        }
//    }
//
//    for (int i = 1; i < node_number; i ++){
//        if (nodes[i] != NONE) {
//
//        }
//    }
//
//    int node_id = 0;
//    for (auto & node_neighbors: edges){
//        // check if all neighbors of node are opposite part
//        int neighbor_id = 0;
//        for (auto neighbor: node_neighbors){
//            if (neighbor && (bipartite_nodes[node_id] == bipartite_nodes[neighbor_id]) && (node_id != neighbor_id))
//                return false;
//            neighbor_id ++;
//        }
//        node_id ++;
//    }
//
//    return true;
//}

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
//
//void Graph::calculate_aux(vector<vector<bool>> subgraph_edges, vector<BIPART> bipartite_nodes, int &maximum,
//                          vector<vector<vector<bool>>> &solutions) {
//    // stop condition
//    if (bipartite_nodes.size() == this->n){
//        int current_weight = get_edges_weight(subgraph_edges);
//        // save best solutions
//        if (current_weight >= maximum){
//            maximum = current_weight;
//            vector<vector<bool>> solution = subgraph_edges;
//            if (current_weight > maximum)
//                solutions.clear();
//            solutions.push_back(solution);
//        }
//        return;
//    }
//
//    // go through neighbors
//    for (int i = 0; i < subgraph_edges[bipartite_nodes.size() - 1].size(); i ++){
//        if (subgraph_edges[bipartite_nodes.size() - 1][i]){
//            // assign to part one
//            vector<BIPART> b_one = bipartite_nodes;
//            b_one.push_back(BIPART::ONE);
//            calculate_aux(, b_one, maximum, solutions);
//
//            // assign to part two
//        }
//    }
//
//}

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

pair<int, int> Graph::get_edge_vertices (int edge_idx) {

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
    bool predefined_vertices = false;
    if ((bipartite_nodes[f] == ONE && bipartite_nodes[s] == TWO) ||
            (bipartite_nodes[f] == TWO && bipartite_nodes[s] == ONE)) {
        predefined_vertices = true;
        calculate_aux(new_subgraph_edges_with, new_nodes, maximum, solutions, ++node_depth,
                      weight + this->edges_weight[f][s]);
    }
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
    if (!predefined_vertices) {
        vector<bool> new_subgraph_edges_without = subgraph_edges;
        new_subgraph_edges_without.push_back(false);
        calculate_aux(new_subgraph_edges_without, bipartite_nodes, maximum, solutions, ++node_depth, weight);
    }
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
//    cout << "is coherent? " << is_coherent(solutions[0]) << endl;

    recursion_count = 0;

    return pair(maximum, (int) solutions.size());
}
