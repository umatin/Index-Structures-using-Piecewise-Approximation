#pragma once

#include <algorithm>
#include <chrono>
#include <iostream>
#include <functional>
#include <fstream>
#include <thread>
#include <vector>
#include <cassert>
#include <unordered_map>


// Loads values from binary file into vector.
template<typename T>
static std::vector<T> load_data(const std::string& filename) {
    std::vector<T> data;
    
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cout << "unable to open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    // Read size.
    uint64_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    data.resize(size);
    
    // Read values.
    in.read(reinterpret_cast<char*>(data.data()), size*sizeof(T));
    in.close();

    return data;
}
