//
// Created by tomas on 18.02.22.
//

#ifndef SEQUENTIAL_GRAPH_H
#define SEQUENTIAL_GRAPH_H

#include <vector>
#include <string>
#include <memory>

using namespace std;

enum BIPART {
    ONE, TWO, NONE
};

class Graph {
private:
    /// number representing number of nodes - 50 > n >= 10
    int n;

    /// number represents average number of node's degree - n/2 >= k >= 3
    int k;

    /// number of edges
    int e;

    /// filepath to
    string filepath;

    /// weighted edges_weight between interval <80,120>
    vector < vector<int> > edges_weight;

    /// reduced representation of edges
    vector<int> edges_weight_redux;

    /// incident vertices with edge
    vector<pair<int, int>> edges_vertices;

    int num_threads;

public:

    explicit Graph(const string & filepath, int num_threads);

    void calculate_aux(unique_ptr<vector<bool>> subgraph_edges, unique_ptr<vector<BIPART>> bipartite_nodes, unique_ptr<int> node_depth, unique_ptr<int> weight);

    void calculate ();

    static bool is_coherent(const vector< vector<bool> >& edges);

    static void is_coherent_aux(const vector<vector<bool>> &edges, int node_idx, vector<bool> &visited);

    vector<vector<bool>> translate_vector_to_edges(const vector<bool> &edges);

    int get_edges_weight(const vector< vector<bool> >& subgraph_edges);

    bool has_potential(const vector<bool> &subgraph_edges, int maximum, int weight);

    pair<int, int> get_edge_vertices(int edge_idx);

    vector<bool> translate_edges_to_vector();

    static bool is_bipartite(vector<BIPART> bipartite_nodes, const vector< vector<bool> >& edges);

};


#endif //SEQUENTIAL_GRAPH_H
