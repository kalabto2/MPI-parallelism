//
// Created by tomas on 18.02.22.
//

#include "Graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>

#define CSV

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

#ifndef CSV
    cout << this->filepath << " -- Solution: " << maximum << " (" << solutions.size() << ") "
        << "(" << duration << " s) (recursion - " << recursion_count << " )" << endl;
#endif
#ifdef CSV
    // output to csv
    string sep = ";";
    // filepath, n, k, solution, # solutions, duration in seconds, running threads, running processes
    cout << this->filepath << sep << this->n << sep << this->e << sep << maximum << sep << solutions.size() <<
    sep << duration << sep << 1 << sep << 1 << endl;
#endif
    recursion_count = 0;

    return pair(maximum, (int) solutions.size());
}
