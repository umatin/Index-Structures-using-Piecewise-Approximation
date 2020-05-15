#pragma once

#include "convex_hull_tree.h"
#include <vector>
#include <iostream>


template<class K>
class segment{
public:
    ch_tree<K> ch;
    double a = 0;
    double b = 0;

    segment(){

    }

    ~segment(){
        
    }

    long predict(K x){
        return (long) ((a*(double)x)+b);
    }

    void insert(K x, long y){
        ch.insert(x,y);
    }

    void remove(K k){
        ch.remove(k);
    }

    void print(){
        ch.printInline();
    }
    
    void printshort(){
        DBGP(Lower)
        std::cout << ch.printshort(ch.rootLower, 0);
        DBGP(Upper)
        std::cout << ch.printshort(ch.rootUpper, 0);
    }

    void split(segment * other, long x){
        ch.split(&(other->ch), x);
    }

    std::pair<tree_node<K> *, tree_node<K>*> refresh_spline(double * width){
        auto ret = ch.get_min_width(width, &a, &b);
        return ret;
    }

    void reset(){
        ch.reset();
    }

    long get_Size(){
        return ch.getSize();
    }
};