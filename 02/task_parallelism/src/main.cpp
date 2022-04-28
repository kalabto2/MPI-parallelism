//
// Created by tomas on 18.02.22.
//

#include <memory>
#include <iostream>
#include "Instance_handler.h"
#include "Graph.h"

using namespace std;

const string TEST_DIRECTORY_PATH = "../graf_mbp/";

int main(int argc, char ** argv) {
    if (argc == 1) {
        cout << "bad number of arguments" << endl;
        return 1;
    }
    int num_threads = stoi(argv[1]);
    if (argc == 2) {
        unique_ptr <Instance_handler> handler(new Instance_handler(TEST_DIRECTORY_PATH, num_threads));
        handler->test_all();
    }
    else if (argc == 3){
        string fp = string(argv[2]);
        unique_ptr<Graph> graph (new Graph(fp, num_threads));
        graph->calculate();
    }


    return 0;
}