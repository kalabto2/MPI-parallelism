//
// Created by tomas on 18.02.22.
//

#include "Instance_handler.h"
#include "Graph.h"

#include <memory>
#include <filesystem>
#include <iostream>

#define CSV

Instance_handler::Instance_handler(const string & dir_path, int num_threads) {
    for (const auto & entry : filesystem::directory_iterator(dir_path)){
        string file_path = entry.path();
        file_paths.push_back(file_path);
    }
    this->num_threads = num_threads;

    correct_solutions["graf_10_3.txt"] = "1300 (2) (0,0 s)\t(recursion - 1K)";
    correct_solutions["graf_10_5.txt"] = "1885 (1) (0,0 s)\t(recursion - 12K)";
    correct_solutions["graf_10_6.txt"] = "2000 (2) (0,02 s)\t(recursion - 374K)";
    correct_solutions["graf_10_7.txt"] = "2348 (1) (0,01 s)\t(recursion - 360K)";
    correct_solutions["graf_12_3.txt"] = "1422 (1) (0,0 s)\t(recursion - 1K)";
    correct_solutions["graf_12_5.txt"] = "2219 (1) (0,0 s)\t(recursion - 171K)";
    correct_solutions["graf_12_6.txt"] = "2533 (1) (0,0 s)\t(recursion - 401K)";
    correct_solutions["graf_12_9.txt"] = "3437 (1) (0,47 s)\t(recursion - 49M)";
    correct_solutions["graf_15_4.txt"] = "2547 (1) (0,0 s)\t(recursion - 502K)";
    correct_solutions["graf_15_5.txt"] = "2892 (1) (0,03 s)\t(recursion - 1.3M)";
    correct_solutions["graf_15_6.txt"] = "3353 (1) (0,05 s)\t(recursion - 4M)";
    correct_solutions["graf_15_8.txt"] = "3984 (1) (1,8 s)\t(recursion - 160M)";
}

void Instance_handler::test_all() {
    for (auto & fp: file_paths){
        unique_ptr<Graph> graph (new Graph(fp, this->num_threads));
        graph->calculate();

#ifndef CSV
        string name = fp.substr(fp.size() - 13, 13);
        std::cout << "\t\t\t    Correct: " << this->correct_solutions[name] << std::endl;
#endif
    }
}
