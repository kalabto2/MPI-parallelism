//
// Created by tomas on 18.02.22.
//

#include <memory>
#include "Instance_handler.h"

using namespace std;

const string TEST_DIRECTORY_PATH = "../graf_mbp/";

int main(int argc, char ** argv) {
    unique_ptr <Instance_handler> handler(new Instance_handler(TEST_DIRECTORY_PATH));
    handler->test_all();

    return 0;
}