#pragma once

#include "convex_hull_tree.h"
#include "segment.h"
#include <vector>
#include <iostream>

template<class K, class V>
class key_val{
public:
    K key;
    V val;
};

template<class K, class V>
class gapped_array{
public:

    union data_ptr {
        key_val<K,V> * store = nullptr;
        key_val<K,gapped_array<K,V> *> * node;
    };
    data_ptr data;
    gapped_array<K, V> * next_ga = nullptr;
    gapped_array<K, V> * parent = nullptr;
    segment<K> spl;
    
    bool isNode = true;
    long s = 0;
    long container_size = 0;
    long base_size;
    double max_epsilon;


    //#define EXISTS(i) (i>0 ? data[i].key != data[i-1].key : true)
    
    gapped_array(long base_size, bool isNode, double max_epsilon){
        data.node = nullptr;
        this->isNode = isNode;
        this->base_size = base_size;
        this->max_epsilon = max_epsilon;

    }

    ~gapped_array(){

    }

    bool get_isNode(){
        return isNode;
    }

    template<class value_type>
    inline bool exists(long i, key_val<K,value_type> * d){
        return (i>0) ? d[i].key != d[i-1].key : true;
    }
    
    inline K get_min_key(){
        if(isNode){
            return data.node[0].key;
        }
        else{
            return data.store[0].key;
        }
    }

    inline void setParent(gapped_array<K, V> * p){
        parent = p;
    }

    void init(){
        this->container_size = base_size;
        if(isNode){
            data.node = new key_val<K, gapped_array<K,V>*>[base_size];
            for(int i = 0; i < base_size; i++){
                data.node[i].key = 0;
                data.node[i].val = nullptr;
            }
        }
        else{
            data.store = new key_val<K, V>[base_size];
            for(int i = 0; i < base_size; i++){
                data.store[i].key = 0;
            }
        }
    }
     
    template<class value_type>
    long next(long pos, key_val<K,value_type> * d){
        K k = d[pos].key;

        while(k == d[pos].key && pos < container_size){
            pos++;
        }

        return pos;
    }
    
    template<class value_type>
    long prev(long pos, key_val<K,value_type> * d){
        K k = d[pos].key;

        do{
            pos--;
        }while(!exists<value_type>(pos, d) && pos >= 0);

        return pos;
    }

    void print(){
        if(isNode){
            DBGP(isNode)
            for (int i = 0; i < container_size; i++){
                if(exists<gapped_array<K,V> * >(i, data.node)){
                    std::cout << data.node[i].key << "\t";
                }
                else{
                    std::cout << "\033[0;31m" << data.node[i].key << "\033[0m\t";
                }
                if((i+1) % 4 == 0){
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
        else{
            DBGP(is NOT Node)
            for (int i = 0; i < container_size; i++){
                if(exists<V>(i, data.store)){
                    std::cout << data.store[i].key << "\t";
                }
                else{
                    std::cout << "\033[0;31m" << data.store[i].key << "\033[0m\t";
                }
                if((i+1) % 4 == 0){
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
        DBGPA(this->s)
    }

    void print_data(){
        DBGPB(this->s, this->container_size)
        if(isNode){
            DBGP(isNode)
            for (int i = 0; i < container_size; i++){
                if(exists<gapped_array<K,V> * >(i, data.node)){
                    std::cout << data.node[i].key << ":" << data.node[i].val << "\t";
                }
                else{
                    std::cout << "\033[0;31m" << data.node[i].key << ":" << data.node[i].val << "\033[0m\t";
                }
                if((i+1) % 4 == 0){
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
        else{
            DBGP(is NOT Node)
            for (int i = 0; i < container_size; i++){
                if(exists<V>(i, data.store)){
                    std::cout << data.store[i].key << ":" << data.store[i].val << "\t";
                }
                else{
                    std::cout << "\033[0;31m" << data.store[i].key << ":" << data.store[i].val << "\033[0m\t";
                }
                if((i+1) % 4 == 0){
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
        
    }

    inline gapped_array<K,V> * get(key_val<K,V> ** ret_kv, K key){
        long prediction = spl.predict(key);
        if(isNode){
            //long pos = get_pos_leq<gapped_array<K,V> *>(key, prediction, data.node);
            long pos = get_pos_leq_exponentialSearch<gapped_array<K,V> *>(key, prediction, data.node);
            return data.node[pos].val;
        }
        else{
            //long pos = get_pos_leq<V>(key, prediction, data.store);
            long pos = get_pos_leq_exponentialSearch<V>(key, prediction, data.store);
            *ret_kv = &data.store[pos];
            return this;
        }
    }

    template<class value_type>
    inline long exponentialSearch(K x, int off, key_val<K,value_type> * d) 
    {
        if (d[off].key == x){
            return off; 
        }
        long i = 1; 
        while (off+i < container_size ? d[off+i].key <= x : false){
            i = i*2; 
        }
        return binarySearch(x, off+(i/2), std::min(off+i+1, container_size), d); 
    }

    template<class value_type>
    inline long exponentialSearchBW(K x, int off, key_val<K,value_type> * d) 
    {
        if (d[off].key == x){
            return off+1; 
        }
        long i = 1; 
        while (off-i > 0 ? d[off-i].key >= x : false){
            i = i*2; 
        }
        return binarySearch(x, std::max(off-i, (long)0), off-(i/2)+1, d); 
    } 
    

    template<class value_type>
    static inline long binarySearch(K x, long fr, long to, key_val<K,value_type> * d){
        long loc;
        long step;
        long count = to - fr;

        while (count > 0) {
            loc = fr; 
            step = count / 2; 
            loc += step;
            if(!(x < d[loc].key)) {
                ++loc;
                fr = loc;
                count -= step + 1;
            }
            else {
                count = step;
            }
        }
        return fr;
    }

    template<class value_type>
    inline long get_pos_leq_exponentialSearch(K x, long prediction, key_val<K,value_type> * d){
        //std::cout << x << std::endl;
        //print_data();
        long i = prediction;
        i = std::max(i,0l);
        i = std::min(i,container_size-1);
        
        if(d[i].key < x){
            long ret = exponentialSearch(x, i, d);
            if(ret >= container_size){
                ret--;
            }
            while(ret > 0 && d[ret].key > x){
                ret--;
            }
            return ret;
        }
        else if(d[i].key > x){
            long ret = exponentialSearchBW(x, i, d);
            if(ret >= container_size){
                ret--;
            }
            while(ret > 0 && d[ret].key > x){
                ret--;
            }
            return ret;
        }
        else{
            return i;
        }
    }

    template<class value_type>
    long get_pos_leq(K x, long prediction, key_val<K,value_type> * d){
        long i = prediction;
        i = std::max(i,0l);
        i = std::min(i,container_size-1);
        
        if(d[i].key < x){
            while(i < container_size){
                if(d[i].key > x){
                    break;
                }
                else{
                    i++;
                }
            }
            i--;

            //assert(d[i].key <= x || i == 0);
            //for(int j = i; j < container_size; j++){
            //    assert(d[i].key == d[j].key || x < d[j].key);
            //}
            return i;
        }
        else if(d[i].key > x){
            while(i > 0){  // this is intentional
                if(d[i].key <= x){
                    break;
                }
                else{
                    i--;
                }
            }
            //assert(d[i].key <= x || i == 0);
            //for(int j = i; j < container_size; j++){
            //    assert(d[i].key == d[j].key || x < d[j].key);
            //}
            return i;
        }
        else{
            return i;
        }
    }
    
    template<class value_type>
    long get_pos(K x, long prediction, key_val<K,value_type> * d){
        long i = prediction;
        i = std::max(i,0l);
        i = std::min(i,container_size-1);
        
        if(d[i].key < x){
            while(i < container_size){
                if(d[i].key >= x){
                    break;
                }
                else{
                    i++;
                }
            }
            if(!exists<value_type>(i-1,d) && d[i].key > x && i != 0){
                return i-1;
            }
            return i;
        }
        else if(d[i].key > x){
            while(i > 0){  // this is intentional
                if(d[i-1].key < x){
                    break;
                }
                else{
                    i--;
                }
            }
            if(!exists<value_type>(i-1,d) && i != 0){
                return i-1;
            }
            return i;
        }
        else{
            return i;
        }
    }

    gapped_array<K,V> * insert_node(K x, gapped_array<K,V> * v){
        if(s == 0){
            init();
        }
        if(s*2 >= container_size){
            spread<gapped_array<K,V> *>();
        }
        insert_base<gapped_array<K,V> *>(x,v,data.node);
        
        //print_data();
        double width;
        auto tree_node_pair = spl.refresh_spline(&width);
        if(width > max_epsilon){
            return split<gapped_array<K,V> *>();
        }
        return nullptr;
    }

    gapped_array<K,V> * insert(K x, V v){
        if(s == 0){
            init();
        }
        if(s*2 >= container_size){
            spread<V>();
        }
        insert_base<V>(x,v,data.store);

        double width;
        auto tree_node_pair = spl.refresh_spline(&width);
        if(width > max_epsilon){
            return split<K>();
        }
        return nullptr;
    }

    void replace_min(K old_min, K new_min){
        assert(new_min < old_min);
        assert(isNode);

        //long pos = get_pos_leq<gapped_array<K,V> *>(old_min, spl.predict(old_min), data.node);
        long pos = get_pos_leq_exponentialSearch<gapped_array<K,V> *>(old_min, spl.predict(old_min), data.node);

        for (int i = pos; i < container_size; i++){
            if(data.node[i].key == old_min){
                data.node[i].key = new_min;
            }
            else{
                break;
            }
        }
        for (int i = pos-1; i >= 0; i--){
            if(data.node[i].key == old_min){
                data.node[i].key = new_min;
            }
            else{
                break;
            }
        }
    }

    //return value is new root or nullptr
    template<class value_type>
    void insert_base(K x, value_type v, key_val<K,value_type> * d){

        if(get_min_key() > x){
            if(parent != nullptr){
                parent->replace_min(get_min_key(), x);
            }
        }

        int i = 0;
        //the first element should always be inserted at 0
        if(s != 0){
            long prediction = spl.predict(x);
            i = get_pos<value_type>(x, prediction, d);

            if(i < container_size){
                if(d[i].key == x){
                    //assert(exists<value_type>(i,d));
                    //TODO change val
                    return ;//already inside
                }
            }

            if (i > 0){
                if(!exists<value_type>(i-1, d)){
                    i -= 1;
                }else{
                    if(i < container_size){
                        if(!shiftRight<value_type>(i, d)){
                            shiftLeft<value_type>(i-1, d);
                            i-=1;
                        }else{
                            //all ok
                        }
                    }else{
                        shiftLeft<value_type>(i-1, d);
                        i-=1;
                    }
                }
            }else{
                shiftRight<value_type>(i, d);
            }
        }
        d[i].key = x;
        d[i].val = v;
        s++;
        spl.insert(x,i);

        if(i > 0){
            K pr_val = d[i-1].key;
            i++;
            while(i < container_size ? d[i].key == pr_val : false){
                d[i].key = x;
                d[i].val = v;
                i++;
            }
        }
        if(s == 1){
            while(i < container_size){
                d[i].key = x;
                d[i].val = v;
                i++;
            }
        }
    }

    template<class value_type>
    void spread(){
        if(s>50){
            //DBGP(spread)
            //DBGPB(s, container_size)
        }
        spl.reset();
        if(isNode){
            create<gapped_array<K,V> *>(2*container_size, container_size, data.node, &data.node);
        }
        else{
            create<V>(2*container_size, container_size, data.store, &data.store);
        }
    }

    //do copy and redistribute favorably
    template<class value_type>
    void create(long target_len, long source_len, key_val<K,value_type> * prev_store, key_val<K,value_type> ** target_store){
        spl.reset();

        *target_store = new key_val<K,value_type>[target_len];
        container_size = target_len;

        double dx = target_len;
        double dy = prev_store[source_len-1].key - prev_store[0].key;
        double a = dx/dy;
        double b = 0 - ( a * prev_store[0].key);

        (*target_store)[0] = prev_store[0];

        long prev_loc = 1;
        key_val<K, value_type> prev_val = prev_store[0];
        for(long i = 1; i < source_len; i++){
            if(exists<value_type>(i,prev_store)){
                int p = (long) ((a*(double)prev_store[i].key)+b);

                if(prev_loc > target_len-1){
                    prev_loc = target_len-1;
                    
                    //Shift left basically
                    long pos = prev_loc;
                    while(exists<value_type>(pos,(*target_store)) && pos >= 0){
                        pos--;
                    }
                    for(; pos < prev_loc; pos++){
                        (*target_store)[pos] = (*target_store)[pos+1];
                    }
                }
                else{
                    if(p > prev_loc){
                        while(prev_loc < p && prev_loc < target_len-1){
                            (*target_store)[prev_loc] = prev_val;
                            prev_loc++;
                        }
                    }
                }
                prev_val = prev_store[i];
                (*target_store)[prev_loc] = prev_store[i];
                prev_loc++;
            }
        }

        while(prev_loc < target_len){
            (*target_store)[prev_loc] = prev_val;
            prev_loc++;
        }

        for(int i = 0; i < target_len; i++){
            if(exists<value_type>(i,(*target_store))){
                spl.insert((*target_store)[i].key, i);
            }
        }

        double w;
        this->spl.refresh_spline(&w);
    }  

    //split into [0, pos) and [pos, container_size)
    template<class value_type>
    gapped_array<K,V> * split(){
        //print_data();
        //DBGP(split)
        //DBGPB(s, container_size)

        gapped_array<K,V> * new_ga = new gapped_array<K,V>(base_size, isNode, max_epsilon);

        data_ptr data_old = data;

        if(isNode){
            long pos = 0;
            long count = 0;
            while (count < s/2 || !exists<gapped_array<K,V> *>(pos, data.node)){
                if(exists<gapped_array<K,V> *>(pos, data.node)){
                    count++;
                }
                pos++;
            }
            if(pos < container_size){
                long new_len = container_size - pos;
                long target_len = container_size/2;
                
                new_ga->create<gapped_array<K,V>*>(target_len, new_len, &data.node[pos], &(new_ga->data.node));
                this->create<gapped_array<K,V>*>(target_len, pos, data.node, &data.node);
                delete data_old.node;
            }
            else{
                assert(false);
                return nullptr;
            }
        }
        else{
            long pos = 0;
            long count = 0;

            while (count < s/2 || !exists<V>(pos, data.store)){
                if(exists<V>(pos, data.store)){
                    count++;
                }
                pos++;
            }

            if(pos < container_size){
                long new_len = container_size - pos;
                long target_len = container_size/2;

                new_ga->create<V>(target_len, new_len, &data.store[pos], &(new_ga->data.store));
                
                this->create<V>(target_len, pos, data.store, &data.store);
                delete data_old.store;
            }
            else{
                assert(false);
                return nullptr;
            }
        }


        //Singly linked list insert
        new_ga->next_ga = this->next_ga;
        this->next_ga = new_ga;
        new_ga->parent = this->parent;

        s = 0;
        //new_ga->print_data();
        update_parent_and_size();
        new_ga->update_parent_and_size();

        if(parent != nullptr){
            assert(parent->isNode);
            return parent->insert_node(new_ga->get_min_key(), new_ga);
        }
        else{
            gapped_array<K,V> * new_root = new gapped_array<K,V> (base_size, true, max_epsilon);
            this->setParent(new_root);
            new_ga->setParent(new_root);
            new_root->insert_node(get_min_key(), this);
            new_root->insert_node(new_ga->get_min_key(), new_ga);
            
            return new_root;
        }
    }

    void update_parent_and_size(){
        if(isNode){
            for(int i = 0; i < container_size; i++){
                if(exists<gapped_array<K,V>*>(i,data.node)){
                    data.node[i].val->setParent(this);
                    s++;
                }
            }
        }
        else{
            for(int i = 0; i < container_size; i++){
                if(exists<V>(i,data.store)){
                    s++;
                }
            }
        }
    }
    
    template<class value_type>
    bool shiftLeft(long i, key_val<K,value_type> * d){
        assert(i < container_size && i >= 0);
        long pos = i;

        while(exists<value_type>(pos, d) && pos >= 0){
            pos--;
        }

        if (pos >= 0){
            for(; pos < i; pos++){
                d[pos] = d[pos+1];

                spl.remove(d[pos].key);
                spl.insert(d[pos].key, pos);
            }
            return true;
        }
        else{
            return false;
        }
    }

    template<class value_type>
    bool shiftRight(long i, key_val<K,value_type> * d){
        assert(i < container_size && i >= 0);
        long pos = i;

        while(exists<value_type>(pos, d) && pos < container_size){
            pos++;
        }
        if (pos < container_size){
            
            for(; i < pos; pos--){
                d[pos] = d[pos-1];
                                
                spl.remove(d[pos].key);
                spl.insert(d[pos].key, pos);
            }
            return true;
        }
        else{
            return false;
        }
    }

    long count(){
        assert(!isNode);

        long cnt=0;
        for(int i = 0; i < container_size; i++){
            if(exists<V>(i,data.store)){
                cnt++;
            }
        }
        return cnt;
    }

    gapped_array<K,V>* check_sorted(K * min_k){
        if(!isNode){
            assert(data.store[0].key > *min_k);
            for(int i = 1; i < container_size; i++){
                assert(data.store[i].key >= data.store[i-1].key);
            }
            *min_k = data.node[container_size - 1].key;
        }
        else{
            assert(data.node[0].key > *min_k);
            for(int i = 1; i < container_size; i++){
                assert(data.node[i].key >= data.node[i-1].key);
            }
            *min_k = data.node[container_size - 1].key;
        }
        return next_ga;
    }

    void print_contains(K x){
        if(!isNode){
            for(int i = 1; i < container_size; i++){
                if(data.store[i].key == x){
                    print_data();
                    DBGP(found)
                    return;
                }
            }
        }
    
    }
};