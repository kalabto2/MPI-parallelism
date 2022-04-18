#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <mpi.h>
#include <memory>
#include <queue>


#define LENGTH 10000
#define TAG_WORK 0
#define TAG_FINISHED 1
#define TAG_TERMINATE 2


//#define DEBUG
//#define DEBUG_PACKING

using namespace std;

enum BIPART {
    ONE, TWO, NONE
};

int position = 0;
char buffer[LENGTH];

int recursion_count = 0;
int maximum = -1;
//vector<vector<vector<bool>>> solutions;

class STATE {
public:
    STATE() {}
    STATE(vector<bool> subgraphEdges, vector<BIPART> bipartiteNodes, int nodeDepth, int weight)
            : subgraph_edges(std::move(subgraphEdges)), bipartite_nodes(std::move(bipartiteNodes)), node_depth(nodeDepth), weight(weight) {}

    vector<BIPART> bipartite_nodes;
    vector<bool> subgraph_edges;
    int node_depth;
    int weight;
};

void packMPI (const vector<STATE> & states, const int global_maximum) {
    int position = 0;

    int states_size = (int)states.size();
    MPI_Pack(&states_size, 1, MPI_INT, buffer, LENGTH, &position, MPI_COMM_WORLD);
    MPI_Pack(&global_maximum, 1, MPI_INT, buffer, LENGTH, &position, MPI_COMM_WORLD);

    for (auto state: states){
#ifdef DEBUG_PACKING
        cout << "PACKING: [";
        cout << state.weight << endl;
        cout << state.node_depth << endl;
        for (auto e : state.subgraph_edges)
            cout << e << " ";
        cout << endl;
        for (auto e : state.bipartite_nodes)
            cout << e << " ";
        cout << endl;
#endif
        MPI_Pack(&state.weight, 1, MPI_INT, buffer, LENGTH, &position, MPI_COMM_WORLD);
        MPI_Pack(&state.node_depth, 1, MPI_INT, buffer, LENGTH, &position, MPI_COMM_WORLD);

        unsigned edges_size = state.subgraph_edges.size();
        int * e_formated = new int [edges_size];
        for (int i = 0; i < edges_size; i ++)
            e_formated[i] = state.subgraph_edges[i];

        MPI_Pack(&edges_size, 1, MPI_UNSIGNED, buffer, LENGTH, &position, MPI_COMM_WORLD);
        MPI_Pack(e_formated, (int)edges_size, MPI_INT, buffer,
                 LENGTH, &position, MPI_COMM_WORLD);

        unsigned nodes_size = state.bipartite_nodes.size();
        int * n_formated = new int [nodes_size];
        for (int i = 0; i < nodes_size; i ++)
            n_formated[i] = state.bipartite_nodes[i];

        MPI_Pack(&nodes_size, 1, MPI_UNSIGNED, buffer, LENGTH, &position, MPI_COMM_WORLD);
        MPI_Pack(n_formated, (int)nodes_size, MPI_INT,
                 buffer, LENGTH, &position, MPI_COMM_WORLD);

        delete[] n_formated;
        delete[] e_formated;
    }

}

void unpackMPI(vector<STATE> & states)
{
    position = 0;
    int new_maximum;
    int state_size;
    MPI_Unpack(buffer, LENGTH, &position, &state_size, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, LENGTH, &position, &new_maximum, 1, MPI_INT, MPI_COMM_WORLD);

#pragma omp critical
    maximum = new_maximum;

    for (int j = 0; j < state_size; j ++) {
        STATE state;

        MPI_Unpack(buffer, LENGTH, &position, &state.weight, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(buffer, LENGTH, &position, &state.node_depth, 1, MPI_INT, MPI_COMM_WORLD);

        unsigned edges_size, node_size;
        MPI_Unpack(buffer, LENGTH, &position, &edges_size, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        int *e_formated = new int[edges_size];
        state.subgraph_edges.resize(edges_size);
        MPI_Unpack(buffer, LENGTH, &position, e_formated, (int) edges_size, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < edges_size; i++)
            state.subgraph_edges[i] = e_formated[i];

        MPI_Unpack(buffer, LENGTH, &position, &node_size, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        int *n_formated = new int[node_size];
        state.bipartite_nodes.resize(node_size);
        MPI_Unpack(buffer, LENGTH, &position, n_formated, node_size, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < node_size; i++)
            state.bipartite_nodes[i] = (BIPART) n_formated[i];

        states.push_back(state);

#ifdef DEBUG_PACKING
        cout << "unpacking: " << endl;
        cout << "--------------" << endl;
        cout << ">> e_size: " << edges_size << endl;
        cout << ">> n_size: " << node_size << endl;
        cout << state.weight << endl;
        cout << state.node_depth << endl;

        for (auto e : state.subgraph_edges)
            cout << e << " ";
        cout << endl;
        for (auto e : state.bipartite_nodes)
            cout << e << " ";
        cout << endl;
#endif

        delete[] n_formated;
        delete[] e_formated;
    }
}


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
    vector<STATE> solution;

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

    void get_n_subgraphs(vector<STATE> &output_states, const vector<STATE> &input_states, int n, int depth);
};

// ================================================================================================

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

    int resting_weight = 0;
    for (int i = curr_e; i < this->e; i ++) {
        resting_weight += this->edges_weight_redux[i];
    }
    return weight + resting_weight >= maximum;
}

void Graph::calculate_aux(const vector<bool>& subgraph_edges, const vector<BIPART>& bipartite_nodes,
                          int node_depth, int weight) {
#pragma omp atomic
    recursion_count ++;

    // stop conditions
    if (subgraph_edges.size() >= this->e){
        vector<vector<bool>> translated_edges = this->translate_vector_to_edges(subgraph_edges);
        int current_weight = get_edges_weight(translated_edges);

        if (current_weight >= maximum &&
            is_coherent(translated_edges) &&
            all_of(bipartite_nodes.begin(), bipartite_nodes.end(), [](BIPART i){return i != NONE;})){
            // save maximum
#pragma omp critical
            if (current_weight >= maximum &&
                is_coherent(translated_edges) &&
                all_of(bipartite_nodes.begin(), bipartite_nodes.end(), [](BIPART i){return i != NONE;})) {
                if (current_weight > maximum)
                    this->solution.clear();
                maximum = current_weight;
                this->solution.emplace_back(subgraph_edges, bipartite_nodes, node_depth, current_weight);
            }
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
        calculate_aux(new_subgraph_edges_with, new_nodes, ++node_depth,
                      weight + this->edges_weight[f][s]);
    }
    else if (bipartite_nodes[f] == NONE &&
             bipartite_nodes[s] != NONE) {
        new_nodes[f] = (bipartite_nodes[s] == ONE ? TWO : ONE);
        calculate_aux(new_subgraph_edges_with, new_nodes, ++ node_depth, weight + this->edges_weight[f][s]);
    }
    else if (bipartite_nodes[f] != NONE &&
             bipartite_nodes[s] == NONE) {
        new_nodes[s] = (bipartite_nodes[f] == ONE ? TWO : ONE);
        calculate_aux(new_subgraph_edges_with, new_nodes, ++ node_depth, weight + this->edges_weight[f][s]);
    }
    else if (bipartite_nodes[f] == NONE &&
             bipartite_nodes[s] == NONE) {
        new_nodes[f] = ONE;
        new_nodes[s] = TWO;
        calculate_aux(new_subgraph_edges_with, new_nodes, ++ node_depth, weight + this->edges_weight[f][s]);
        new_nodes[s] = ONE;
        new_nodes[f] = TWO;
        calculate_aux(new_subgraph_edges_with, new_nodes, ++ node_depth, weight + this->edges_weight[f][s]);
    }

    // 2) dont take edge
    if (!predefined_vertices) {
        vector<bool> new_subgraph_edges_without = subgraph_edges;
        new_subgraph_edges_without.push_back(false);
        calculate_aux(new_subgraph_edges_without, bipartite_nodes, ++node_depth, weight);
    }
}

void Graph::calculate() {
    int proc_rank, num_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Status status;

    // ----------------
    vector<STATE> solutions;
    int solution_weight = -1;

    // ================ MASTER ================
    if (proc_rank == 0) {
        // prepare graph
        vector<bool> subgraph_edges;
        vector<BIPART> bipartite_nodes;
        bipartite_nodes.assign(this->n, BIPART::NONE);
        bipartite_nodes[0] = BIPART::ONE;

        // retrieve n disjunctive sub-graphs
        vector<STATE> queue;
        vector<STATE> init_state;
        init_state.emplace_back(subgraph_edges, bipartite_nodes, 0, 0);
        get_n_subgraphs(queue, init_state, 200, 0);

        ::queue<int> free = {};

        // send to all slaves work to do
        for (int i = 1; i < num_processes; i ++){
            STATE job = queue.back();
            queue.pop_back();
            vector<STATE> states;
            states.push_back(job);

            // pack
            packMPI(states, solution_weight);

            // send
            MPI_Send(buffer, LENGTH, MPI_PACKED, i, TAG_WORK, MPI_COMM_WORLD);
#ifdef DEBUG
            cout << "M: send job" << endl;
#endif
        }

        int num_working_slaves = num_processes - 1;

        while (true)
        {
            if (queue.empty() && free.size() == num_working_slaves)
                break;
            int flag = 0;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

            // read received data
            if (flag && status.MPI_TAG == TAG_FINISHED)
            {
                MPI_Recv(buffer, LENGTH, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                vector<STATE> states;
                // unpack
                unpackMPI(states);

                // send new task if available
                if (!queue.empty()) {
                    STATE job = queue.back();
                    queue.pop_back();
                    vector<STATE> states;
                    states.push_back(job);

                    // pack
                    packMPI(states, maximum);
                    // send
                    MPI_Send(buffer, LENGTH, MPI_PACKED, status.MPI_SOURCE, TAG_WORK, MPI_COMM_WORLD);
#ifdef DEBUG
                    cout << "M: sent job" << endl;
#endif
                } else {
#ifdef DEBUG
                    cout << "M: no job lefts" << endl;
#endif
                    MPI_Send(buffer, 0, MPI_C_BOOL, status.MPI_SOURCE, TAG_TERMINATE, MPI_COMM_WORLD);
                    free.push(status.MPI_SOURCE);
                }

                // check if global maximum and save
                if (solution_weight <= maximum) {
                    if (solution_weight <= maximum)
                        this->solution.clear();
                    for (const auto& s: states) {
                        this->solution.push_back(s);
                    }
                    solution_weight = maximum;
                }

            }
        }
#ifdef DEBUG
        cout << "GLOBAL MAXIMUM: " << solution_weight << endl;
        cout << "M: master ends" << endl;
#endif
    } // ================ SLAVE ================
    else {
        while (true)
        {
//            this->solution.clear();

            MPI_Recv(buffer, LENGTH, MPI_PACKED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#ifdef DEBUG
            cout << "S" << proc_rank << ": recieved " << endl;
#endif
            if (status.MPI_TAG == TAG_TERMINATE)
                break;
            else if (status.MPI_TAG == TAG_WORK)
            {
                vector<STATE> states;
                unpackMPI(states);
                STATE job = states[0];

//                cout << "sss " << states.size() << endl;
//                cout << "-------------------" << endl;
//                cout << states[0].node_depth << " " << states[0].weight << endl;
//                cout << states[0].subgraph_edges.size() << endl;
//                for (auto e : states[0].subgraph_edges)
//                    cout << e << "|";
//                cout << endl;
//                for (auto e : states[0].bipartite_nodes)
//                    cout << e << "|";
//                cout << endl;
//                cout << "-------------------" << endl;

                // retrieve n disjunctive sub-graphs
//                vector<STATE> queue;
//                vector<STATE> init_state;
//                init_state.push_back(job);
//                get_n_subgraphs(queue, init_state, 10, job.node_depth);

                calculate_aux(job.subgraph_edges,
                                  job.bipartite_nodes,
                                  job.node_depth,
                                  job.weight);

//                // compute
//#pragma omp parallel for default(none) shared(queue) schedule(dynamic) num_threads(1)
//                for (auto & starting_state : queue) {
//                    calculate_aux(starting_state.subgraph_edges,
//                                  starting_state.bipartite_nodes,
//                                  starting_state.node_depth,
//                                  starting_state.weight);
//                }

                // returning solutions
//                STATE slave_solution;
//                slave_solution.weight = maximum;
//                slave_solution.node_depth = 0;  // not important
//                slave_solution.bipartite_nodes = ;
//                slave_solution.subgraph_edges = translate(solutions[0]);

#ifdef DEBUG
//                cout << "S" << proc_rank << ": n_states: " << queue.size() << endl;
                cout << "S" << proc_rank << ": calculated maximum: " << maximum << endl;
#endif

                packMPI(this->solution, maximum);

                // send to master
                MPI_Send(buffer, LENGTH, MPI_PACKED, 0, TAG_FINISHED, MPI_COMM_WORLD);
            }
        }
#ifdef DEBUG
        cout << "S" << proc_rank << ": ends" << endl;
#endif
    }
}

void Graph::get_n_subgraphs(vector<STATE> &output_states, const vector<STATE> &input_states, int n, int depth) {
    // end if done
    if (depth > this->n)
        return;

    vector<STATE> states;
    for (const auto& input_state : input_states){
        int i = input_state.subgraph_edges.size();
        int f = this->edges_vertices[i].first;
        int s = this->edges_vertices[i].second;

        // 1) take edge
        vector<bool> new_subgraph_edges_with = input_state.subgraph_edges;
        new_subgraph_edges_with.push_back(true);
        vector<BIPART> new_nodes = input_state.bipartite_nodes;
        bool predefined_vertices = false;
        if ((input_state.bipartite_nodes[f] == ONE && input_state.bipartite_nodes[s] == TWO) ||
            (input_state.bipartite_nodes[f] == TWO && input_state.bipartite_nodes[s] == ONE)) {
            predefined_vertices = true;
            states.emplace_back(new_subgraph_edges_with, new_nodes, input_state.node_depth + 1,
                                input_state.weight + this->edges_weight[f][s]);
        }
        else if (input_state.bipartite_nodes[f] == NONE &&
                 input_state.bipartite_nodes[s] != NONE) {
            new_nodes[f] = (input_state.bipartite_nodes[s] == ONE ? TWO : ONE);
            states.emplace_back(new_subgraph_edges_with, new_nodes, input_state.node_depth + 1, input_state.weight + this->edges_weight[f][s]);
        }
        else if (input_state.bipartite_nodes[f] != NONE &&
                 input_state.bipartite_nodes[s] == NONE) {
            new_nodes[s] = (input_state.bipartite_nodes[f] == ONE ? TWO : ONE);
            states.emplace_back(new_subgraph_edges_with, new_nodes, input_state.node_depth + 1, input_state.weight + this->edges_weight[f][s]);
        }
        else if (input_state.bipartite_nodes[f] == NONE &&
                 input_state.bipartite_nodes[s] == NONE) {
            new_nodes[f] = ONE;
            new_nodes[s] = TWO;
            states.emplace_back(new_subgraph_edges_with, new_nodes, input_state.node_depth + 1, input_state.weight + this->edges_weight[f][s]);
            new_nodes[s] = ONE;
            new_nodes[f] = TWO;
            states.emplace_back(new_subgraph_edges_with, new_nodes, input_state.node_depth + 1, input_state.weight + this->edges_weight[f][s]);
        }

        // 2) dont take edge
        if (!predefined_vertices) {
            vector<bool> new_subgraph_edges_without = input_state.subgraph_edges;
            new_subgraph_edges_without.push_back(false);
            states.emplace_back(new_subgraph_edges_without, input_state.bipartite_nodes, input_state.node_depth + 1, input_state.weight);
        }
    }

    // end if enough states already, else continue
    if (states.size() >= n)
        output_states = move(states);
    else
        get_n_subgraphs(output_states, states, n, ++ depth);
}


int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    int proc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    string fp = string(argv[1]);
    unique_ptr<Graph> graph (new Graph(fp));

    // start time
    clock_t start;
    start = clock();
    auto t1 = std::chrono::high_resolution_clock::now();

    graph->calculate();

    // stop time - calculate duration in seconds
    auto t2 = std::chrono::high_resolution_clock::now();
    double duration = (double)( clock() - start ) / (double) CLOCKS_PER_SEC;
    double duration_realtime = (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000.0;

    // print output
    if (proc_rank == 0){
        cout << "------------------------" << endl;
        cout << "FILEPATH: " << fp << endl;
        cout << "DURATION: " << duration << " / " << duration_realtime << endl;
        cout << "SOLUTION: " << graph->solution[0].weight << " (" << graph->solution.size() << ")" << endl;
        cout << "------------------------" << endl;
    }

    MPI_Finalize();
    return 0;
}
