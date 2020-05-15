
#pragma once
#include <vector>
#include <iostream>

#include "convex_hull_tree.h"
#include "gapped_array.h"
#include "segment.h"

template<class K, class V>
class dynamicPWAAindex{
public:
    gapped_array<K,V> * root = nullptr;
    gapped_array<K,V> * store_ll_root = nullptr;
    long size = 0;
    long height = 1;
    long base_size = 0;
    double epsilon = (double)0;

    dynamicPWAAindex(long base_size, double epsilon){
        this->base_size = base_size;
        this->epsilon = epsilon;
        root = new gapped_array<K,V>(base_size, false, (double)epsilon);
        store_ll_root = root;
    }

    template<bool PRINT = false>
    V get(K k){
        gapped_array<K,V> * node = root;
        key_val<K,V> * ret;
        while(node->get_isNode()){
            if(PRINT){
                node->print_data();
            }
            node = node->get(&ret, k);
        }
        node = node->get(&ret, k);
        if(PRINT){
            node->print_data();
        }
        return ret->val;
    }

    void insert(K k, V v){
        gapped_array<K,V> * node = root;
        key_val<K,V> * ret;
        while(node->get_isNode()){
            //node->print_data();
            node = node->get(&ret, k);
        }
        //DBGP(store:)
        //node->print_data();
        node = node->insert(k,v);
        if (node != nullptr){
            root = node;
            this->height++;
        }
    }
    long get_height(){
        return this->height;
    }
    long count(){
        gapped_array<K,V> * node = store_ll_root;
        long cnt=0;
        
        while(node!=nullptr){
            cnt += node->count();
            node = node->next_ga;
        }

        return cnt;
    }

    long count_base_segments(){
        gapped_array<K,V> * node = store_ll_root;
        long cnt=0;
        
        while(node!=nullptr){
            cnt++;
            node = node->next_ga;
        }

        return cnt;
    }

    void check_sorted(){
        gapped_array<K,V> * node = store_ll_root;
        K min_k = -1;
        while(node!=nullptr){
            node = node->check_sorted(&min_k);
        }
    }

    void print_contains(K x){
        gapped_array<K,V> * node = store_ll_root;

        while(node!=nullptr){
            node->print_contains(x);
            node = node->next_ga;
        }
    }

};