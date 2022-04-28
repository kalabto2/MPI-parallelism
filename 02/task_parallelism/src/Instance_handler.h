//
// Created by tomas on 18.02.22.
//

#ifndef SEQUENTIAL_INSTANCE_HANDLER_H
#define SEQUENTIAL_INSTANCE_HANDLER_H

#include <string>
#include <vector>
#include <map>

using namespace std;

class Instance_handler {
private:
    vector<string> file_paths;
    map<string, string> correct_solutions;
    int num_threads;

public:
    Instance_handler(const string & dir_path, int num_threads);

    void test_all ();
};


#endif //SEQUENTIAL_INSTANCE_HANDLER_H
