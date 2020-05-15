#pragma once
#include <vector>
#include <string>
#include <ios>
#include <fstream>
#include <iostream>
#include <stdlib.h>

template<class KeyType>
struct KeyValue {
  KeyType key;
  u_int64_t value;
} __attribute__((packed));

#define RAND(min, max) min + rand() % (( max + 1 ) - min)

//Ease of use Debug print
// TODO remove when #define NDEBUG
#ifdef NDEBUG
    #define DBGP(a) ;
    #define DBGPA(a) ;
    #define DBGPB(a,b) ;
    #define DBGPC(a,b,c) ;
#else
    #define DBGP(a) std::cout << #a << "\n";
    #define DBGPA(a) std::cout << #a << "=" << a << "\n";
    #define DBGPB(a,b) std::cout << #a << "=" << a << " " << #b << "=" << b << "\n";
    #define DBGPC(a,b,c) std::cout << #a << "=" << a << " " << #b << "=" << b << " " << #c << "=" << c << "\n";

#endif
//Avoid double evaluation
#define myMIN(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
#define myMAX(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })

// Loads values from binary file into vector.
template<typename T>
static std::vector<T> load_data(const std::string& filename, bool print = true) {
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

/*
z component of the cross product
or scalar prod of 
(a-o) and (b-0)(flipped by 90 degrees)
Using long is possible too
--->>> y1/x1 > y2/x2 if cross > 0
*/
//no Origin
double cross(double x1, double x1_p, double x2, double x2_p){
    return (x1*x2_p) - (x2*x1_p);
}
//Single origin
double cross(double o, double p_o, double a, double p_a, double b, double p_b){
    return ((a - o) * (p_b - p_o)) - ((p_a - p_o) * (b - o));
}
//Double origin
double cross(double x1_o, double x1_o_p, double x1, double x1_p, double x2_o, double x2_o_p, double x2, double x2_p){
    return ((x1 - x1_o) * (x2_p - x2_o_p)) - ((x1_p - x1_o_p) * (x2 - x2_o));
}