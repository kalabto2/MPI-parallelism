//
// Created by tomas on 18.02.22.
//

#ifndef SEQUENTIAL_GRAPH_H
#define SEQUENTIAL_GRAPH_H

#include <vector>
#include <string>

using namespace std;

enum BIPART {
    ONE, TWO, NONE
};

class STATE {
public:
    STATE(vector<bool> subgraphEdges, vector<BIPART> bipartiteNodes, int nodeDepth, int weight)
            : subgraph_edges(subgraphEdges), bipartite_nodes(bipartiteNodes), node_depth(nodeDepth), weight(weight) {}

    const vector<BIPART> bipartite_nodes;
    const vector<bool> subgraph_edges;
    int node_depth;
    int weight;
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

public:
    explicit Graph(const string & filepath);

    void calculate_aux(const vector<bool>& subgraph_edges, const vector<BIPART>& bipartite_nodes, int node_depth, int weight);

    void calculate ();

    static bool is_coherent(const vector< vector<bool> >& edges);

    static void is_coherent_aux(const vector<vector<bool>> &edges, int node_idx, vector<bool> &visited);

    vector<vector<bool>> translate_vector_to_edges(const vector<bool> &edges);

    int get_edges_weight(const vector< vector<bool> >& subgraph_edges);

    bool has_potential(const vector<bool> &subgraph_edges, int maximum, int weight);

    vector<bool> translate_edges_to_vector();

    static bool is_bipartite(vector<BIPART> bipartite_nodes, const vector< vector<bool> >& edges);

    void get_n_subgraphs(vector<STATE> &output_states, const vector<STATE> &input_states, int n);
};


#endif //SEQUENTIAL_GRAPH_H
