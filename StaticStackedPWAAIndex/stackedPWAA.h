
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <cmath>
#include <cassert>
#include <map>
#include <csignal>
#include <iomanip>

#include <fstream>

//#include "my_util.h"
#include "pwaa_create.h"

template<class K, class V, bool useExponentialSearch = false>
class stackedPWAA{
    public:
    
    std::vector<key_val<K,V>> store;
    std::vector<std::vector<key_val<K,SplineSegment>>> spline_stack;

    double EPSILON;

    stackedPWAA(double e){
        this->EPSILON = e;
    }

    inline long clamp(const long &v, const long &lo, const long &hi){
        return (v < lo) ? lo : ((hi < v) ? hi : v);
    }

    V get(K key){
        long prediction = spline_stack[spline_stack.size()-1][0].val.predict((double)key);
        for(int i = spline_stack.size()-2; i >= 0; i--){
            long pos;
            if(useExponentialSearch){
                pos = get_pos_leq_exponentialSearch(key,prediction,&spline_stack[i][0], spline_stack[i].size());
            }
            else{
                pos = get_pos_leq(key,prediction,&spline_stack[i][0], spline_stack[i].size());
            }
            prediction = spline_stack[i][pos].val.predict((double)key);
        }
        long pos;
        if(useExponentialSearch){
            pos = get_pos_leq_exponentialSearch(key, prediction, &store[0], store.size());
        }
        else{
            pos = get_pos_leq(key, prediction, &store[0], store.size());
        }
        return store[pos].val;
    }
    
    void build(std::vector<K> & d){
        std::vector<key_val<K,V>> data;
        for(int i = 0; i < d.size(); i++){
            data.push_back(key_val<K,V>(d[i],d[i]));
        }
        build_(data);
    }

    void build_(std::vector<key_val<K,V>> & d){
        store = d;
        {
            long fromLoc = 0;
            long toLoc = store.size();
            spline_stack.push_back(std::vector<key_val<K,SplineSegment>>());
            while(fromLoc < toLoc){
                SplineSegment seg;
                K k;
                fromLoc = getSegment<K,V>(fromLoc, toLoc, EPSILON, &k, &seg, &store[0]);
                spline_stack[0].push_back(key_val<K,SplineSegment>(k,seg));
            }
        }


        while(spline_stack[spline_stack.size()-1].size() > 1){            
            spline_stack.push_back(std::vector<key_val<K,SplineSegment>>());
            long fromLoc = 0;
            long toLoc = spline_stack[spline_stack.size()-2].size();
            while(fromLoc < toLoc){
                SplineSegment seg;
                K k;
                fromLoc = getSegment<K,SplineSegment>(fromLoc, toLoc, EPSILON, &k, &seg, &spline_stack[spline_stack.size()-2][0]);
                spline_stack[spline_stack.size()-1].push_back(key_val<K,SplineSegment>(k,seg));
            }
        }
    }

    void print(){
        std::cout << "store size: " << store.size() << std::endl;
        long last_sz = store.size();
        for(int i = 0; i < spline_stack.size(); i++){
            std:: cout << "layer: " << i << " spl size: " << spline_stack[i].size() << std::endl;
            std:: cout << "layer: " << i << " fanout: " << last_sz / spline_stack[i].size() << std::endl;
            last_sz = spline_stack[i].size();
        }

    }
};