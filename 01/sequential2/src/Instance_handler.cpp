//
// Created by tomas on 18.02.22.
//

#include "Instance_handler.h"
#include "Graph.h"

#include <memory>
#include <filesystem>


Instance_handler::Instance_handler(const string & dir_path) {
    for (const auto & entry : filesystem::directory_iterator(dir_path)){
        string file_path = entry.path();
        file_paths.push_back(file_path);
    }
}

void Instance_handler::test_all() {
    for (auto & fp: file_paths){
        unique_ptr<Graph> graph (new Graph(fp));
        graph->calculate();
    }
}
