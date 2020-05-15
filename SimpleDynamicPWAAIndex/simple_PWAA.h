#pragma once

#include "my_util.h"

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



template<class KeyType, bool PMA = false>
class simple_PWAA{
public:
    struct chExtremePoint{
        KeyType x;
        bool exists = true;
    };
    struct SplineSegment {
        long fromY;
        long toY;
        //ax+b
        double a;
        double b;
        double e;
        bool fullSegment;
    };// __attribute__((packed))

    long EPSILON = 64;
    long MAX_ERROR = 64;

    std::vector<KeyValue<KeyType>> store;
    std::vector<bool> exists;
    std::map<KeyType, SplineSegment> spline = {{0,{0,3,0,1,false}}};
    std::map<long, KeyType> spline_segments = {{0,0}};
    
    long levels = 0; //inclusive
    long node_size = 4;
    long container_size = 0;

    //EXPOSED FUNCTIONS

    simple_PWAA(){
        store.resize(node_size);
        exists.resize(node_size);
    }
    
    std::size_t size() const{
        //TODO
        return spline.size() * sizeof(SplineSegment) + store.size() * sizeof(KeyValue<KeyType>);
    }

    void Build(const std::vector<KeyValue<KeyType>>& myData) {
        long segments = (myData.size() / node_size) + ((myData.size() % node_size > 0) ? 1 : 0);
        levels = (long)std::log2(segments) + 1;
        long store_size = (long)pow(2, (levels)) * node_size;
        store.resize(store_size);
        exists.resize(store_size);
        

        assert(store.size() >= myData.size());
        for (long i = 0; i < (long)myData.size(); i++){
            store[i] = myData[i];
            exists[i] = true;
        }
        
        //DBGPC(store_size, segments, levels)
        rebalance(0);
    }

    void insert(KeyValue<KeyType> elem){
        long pos = lookup(elem.key);
        insert_at(elem, pos);

        if (container_size <= 1){
            updateSpline(0, spline.size()-1);
        }

        assert(intact());
        assert(intact_distance());
    }

    KeyValue<KeyType> * get(KeyType x){
        typename std::map<KeyType, SplineSegment>::const_iterator it = spline.upper_bound(x);
        if (it == spline.begin()){
            return nullptr;
        }
        it--;
        SplineSegment spline_seg = it->second;

        DBGPB(spline_seg.a,spline_seg.b)
        long loc = spline_seg.a * (double)x + spline_seg.b;

        if (loc >= (long)store.size()){
            loc = store.size()-1;
        }
        if (loc < 0){
            loc = 0;
        }

        if (store[loc].key == x){
            DBGP(found)
            DBGPB(loc,x)
            return &store[loc];
        }
        else{
            //search to the right
            //TODO <= or < not sure
            for (long i = loc+1; i <= loc+EPSILON+1 && i < (long)store.size(); i++){
                if(exists[i]){
                    if (store[i].key == x){
                        return &store[i];
                    }
                    if (store[i].key > x){
                        break;
                    }
                }
            }
            //search to the left
            for (long i = loc-1; i >= 0 && i > loc-EPSILON-1; i--){
                if(exists[i]){
                    if (store[i].key == x){
                        return &store[i];
                    }
                    if (store[i].key < x){
                        break;
                    }
                }
            }
        }
        return nullptr;
    }
    
    //UTILITY FUNCTIONS
    void checkDistance(long from, long to){
        DBGP(checkDistance)
        DBGPB(from, to)
        auto spl_seg = spline_segments.upper_bound(from);
        assert(spl_seg != spline_segments.begin());
        spl_seg--;
        auto spl = spline.find(spl_seg->second);
        assert(spl != spline.end());
        
        std::vector<std::pair<long,long>> to_refresh;

        long fr = from;
        while (fr <= to){
            long to_value = std::min(to, spl->second.toY);
            for(; fr <= to_value ; fr++){
                double pos = (spl->second.a * store[fr].key) + spl->second.b;
                long addedError = std::abs(fr - (long)pos);
                if(addedError > MAX_ERROR){
                    
                    if (to_refresh.size()>0 ? to_refresh[to_refresh.size()-1].second == spl->second.fromY-1 : false){
                        if(std::next(spl) != spline.end()){
                            to_refresh[to_refresh.size()-1].second = std::next(spl)->second.toY;
                        }
                    }
                    else{
                        if(std::next(spl) != spline.end()){
                            to_refresh.push_back({spl->second.fromY,std::next(spl)->second.toY});
                        }
                    }
                    if(std::next(spl) != spline.end()){
                        spl++;
                        fr = spl->second.toY+1;
                    }else{
                        fr = spl->second.toY+1;
                    }
                }
            }
            spl++;
            spl++;
            //assert(spl->second.fromY == fr);
        }
        for (auto x : to_refresh){
            updateSpline(x.first,x.second);
        }
    }

    //i inclusive
    inline void adjust_spline_forward(long pos, long i){
        //DBGP(adjust_spline_forward)
        //DBGPB(pos, i)

        
        auto spl_seg = spline_segments.upper_bound(pos);
        assert(spl_seg != spline_segments.begin());
        spl_seg--;

        auto spl = spline.find(spl_seg->second);
        assert(spl != spline.end());

        if(spl->second.fromY == pos){
            KeyType tmp_key = store[pos].key;
            spl_seg->second = tmp_key;

            auto nodeHandler = spline.extract(spl);
            nodeHandler.key() = tmp_key;
            auto r = spline.insert(std::move(nodeHandler));
            spl = r.position;
        }

        if(spl->second.toY < i){
            spl->second.toY += 1;
            std::vector<std::pair<long,KeyType>> to_add_segments;
            spl++;
            spl_seg++;
            while(spl != spline.end() && spl->second.toY < i){
                assert(spl->second.fromY == spl_seg->first);
                assert(spl->first == spl_seg->second);
                spl->second.fromY += 1;
                spl->second.toY += 1;

                long tmp_ind = spl_seg->first + 1;
                KeyType tmp_value = spl_seg->second;
                spl_seg = spline_segments.erase(spl_seg);
                to_add_segments.push_back(std::make_pair(tmp_ind,tmp_value));

                spl++;
            }
            if(spl != spline.end() && spl->second.fromY <= i){
                spl->second.fromY += 1;

                long tmp_ind = spl_seg->first + 1;
                KeyType tmp_value = spl_seg->second;
                spl_seg = spline_segments.erase(spl_seg);
                to_add_segments.push_back(std::make_pair(tmp_ind,tmp_value));
            }
            for(auto x : to_add_segments){
                spline_segments.emplace(x.first, x.second);
            }
        }
        assert_print(intact());
    }

    //i inclusive
    inline void adjust_spline_reverse(long pos, long i){
        //DBGP(adjust_spline_reverse)
        //DBGPB(pos, i)


        auto spl_seg = spline_segments.upper_bound(pos);
        assert(spl_seg != spline_segments.begin());
        spl_seg--;

        auto spl = spline.find(spl_seg->second);
        assert(spl != spline.end());

        if(spl != spline.end() && spl->second.toY < i){
            spl->second.toY -= 1;
            
            std::vector<std::pair<long,KeyType>> to_add_segments;
            spl++;
            spl_seg++;

            while(spl->second.toY <= i && spl != spline.end()){
                assert(spl->second.fromY == spl_seg->first);
                assert(spl->first == spl_seg->second);
                spl->second.fromY -= 1;
                if(spl->second.toY < store.size()-1){
                    spl->second.toY -= 1;
                }

                long tmp_ind = spl_seg->first - 1;
                KeyType tmp_value = spl_seg->second;
                spl_seg = spline_segments.erase(spl_seg);
                to_add_segments.push_back(std::make_pair(tmp_ind,tmp_value));

                spl++;
            }

            if(spl != spline.end() && spl->second.fromY == i+1){

                KeyType tmp_key = store[i].key;
                auto tmp_pair = spline.emplace(tmp_key, SplineSegment());
                auto del_spl = spl;
                spl = tmp_pair.first;
                spl->second = del_spl->second;
                spl->second.fromY = i;
                
                spline.erase(del_spl);
                assert(tmp_pair.second);

                long tmp_ind = i;
                KeyType tmp_value = store[i].key;
                spl_seg = spline_segments.erase(spl_seg);
                to_add_segments.push_back(std::make_pair(tmp_ind,tmp_value));

            }
            else if(spl != spline.end() && spl->second.fromY <= i){
                spl->second.fromY -= 1;

                long tmp_ind = spl_seg->first - 1;
                KeyType tmp_value = spl_seg->second;
                spline_segments.erase(spl_seg);
                to_add_segments.push_back(std::make_pair(tmp_ind,tmp_value));
            }
            
            for(auto x : to_add_segments){
                spline_segments.emplace(x.first, x.second);
            }
        }
        assert(intact());
    }

    void vizSlopes(double sl, double off, double err, double err2, double le){
        std::cout << std::setprecision(12) << 0 << " " << off << " " << le << " " << (le*sl) + off << "\n";
        std::cout << std::setprecision(12) << 0 << " " << off + err << " " << le << " " << (le*sl) + off + err << "\n";
        std::cout << std::setprecision(12) << 0 << " " << off + err2 << " " << le << " " << (le*sl) + off + err2 << "\n";
    }

    
    //insert element at chosen location
    void insert_at(KeyValue<KeyType> elem, long pos){
        //Catch invalid calls
        //DBGPB(pos, store.size())
        assert(pos < (long)store.size() && pos >= 0);
        if (store[pos].key == elem.key && exists[pos]) {
            store[pos] = elem;
            return;
        }
        else if(!exists[pos]){
            container_size++;
            store[pos] = elem;
            exists[pos] = true;

            auto spl_seg = spline_segments.upper_bound(pos);
            assert(spl_seg != spline_segments.begin());
            spl_seg--;
            auto spl = spline.find(spl_seg->second);
            assert(spl != spline.end());

            KeyType tmp_key = store[pos].key;
            if (spl_seg->second > tmp_key){
                spl_seg->second = tmp_key;
                spline.emplace(tmp_key, spl->second);
                spline.erase(spl);
            }
            checkDistance(pos, pos);
            rebalance(pos);
            return;
        }
        else{
            container_size++;
            //Find empty slot to the right
            long i = pos;
            while (i < (long)store.size() && exists[i])
            {
                i++;
            }
            if (i < (long)store.size()) {
                //Shift and insert to the right
                
                long write_pos = i;
                exists[write_pos] = true;
                
                while (write_pos > pos) {
                    store[write_pos] = store[write_pos-1];
                    write_pos--;
                }
                store[pos] = elem;

                adjust_spline_forward(pos, i);
                checkDistance(pos, i);

                assert(intact());

                //Potentially 2 Chunks have to be rebalanced
                rebalance(pos);
                rebalance(i);
            }
            else{
                if (store[pos].key > elem.key){
                    pos -= 1;
                }
                i = pos;
                while (i > 0 && exists[i])
                {
                    i--;
                }

                exists[i] = true;
                long write_pos = i;
                while (write_pos < pos){
                    store[write_pos] = store[write_pos+1];
                    write_pos++;
                }
                store[pos] = elem;

                adjust_spline_reverse(i, pos);
                checkDistance(i, pos);

                assert(intact());

                rebalance(pos);
                rebalance(i);
            }
        }
    }

    //finds either an empty spot or the first element greater than x
    long lookup(KeyType x){
        typename std::map<KeyType, SplineSegment>::const_iterator it = spline.upper_bound(x);
        if(it != spline.begin()){
            it--;
        }

        SplineSegment spline_seg = it->second;

        long loc = spline_seg.a * (double)x + spline_seg.b;

        if (loc >= (long)store.size()){
            loc = store.size()-1;
        }
        if (loc < 0){
            loc = 0;
        }

        if (store[loc].key == x && exists[loc]){
            return loc;
        }
        else{
            long i = loc;
            while(i < store.size() ? (!exists[i]) : false){
                i++;
            }
            if (i == store.size()){
                i = loc;
                while(i >= 0 ? !exists[i] : false){
                    i--;
                }
                if (i == -1){
                    return 0;
                }
            }
            assert(exists[i]);
            if(store[i].key > x){
                int i = loc-1;
                while(i >= 0 ? (store[i].key > x || !exists[i]) : false){
                    i--;
                }
                if(store[i].key == x && i >= 0){
                    return i;
                }
                return i+1;
            }
            else if(store[i].key < x){
                int i = loc+1;
                while(i < store.size() ? (store[i].key < x || !exists[i]) : false){
                    i++;
                }

                if(i < store.size()){
                    return i;
                }
                else{
                    //Change
                    return i-1;
                }
            }
            return i;
        }
    }

    

    //PMA UTILITY FUNCTIONS

    //level is distance from root 
    //limits get tighter as they get closer to the root
    inline double lower_limit_density(long level){
        if (levels == 0)
        {
            return 0.25;
        }
        return 0.5 - (0.25 * ((double)level / (double)levels));
    }

    inline double upper_limit_density(long level){
        if (levels == 0)
        {
            return 1.0;
        }
        return 0.75 + (0.25 * ((double)level / (double)levels));
    }

    inline long find_node_ind(long pos){
        return pos / node_size;
    }

    long find_balanced_level(long pos, long * node_ind, long * fr_o, long * to_o){
        *node_ind = find_node_ind(pos);

        //level is distance from root 
        for (long level = levels; level >= 0; level--)
        {
            long dist_from_leaf = levels - level;
            long amt_nodes = pow(2,dist_from_leaf);
            long fr = *node_ind / amt_nodes;//inclusive
            long to = fr + 1;//exclusive
            long amt = 0; //for the to be inserted element
            for (long i = fr*amt_nodes*node_size; i < to*amt_nodes*node_size; i++){
                amt += exists[i] ? 1 : 0;
            }
            double density = (double)amt / ((double)amt_nodes * node_size);
            //DBGPC(density, lower_limit_density(level), upper_limit_density(level))
            //DBGPC(amt_nodes*node_size, amt, density)
            //DBGPB( fr*amt_nodes*node_size, to*amt_nodes*node_size)
            if (density >= lower_limit_density(level) && density <= upper_limit_density(level) && amt < amt_nodes*node_size)
            {
                *fr_o =  fr*amt_nodes*node_size;
                *to_o = (to*amt_nodes*node_size);
                //DBGPB(*fr_o, *to_o)
                return dist_from_leaf;
            }
        }
        // expand has to be called
        return -1;
    }

    //fr inclusive, to exclusive
    long pack(long from, long to){
        //packing stage
        long read_i = from;
        long write_i = from;
        while (read_i < to) {
            if (exists[read_i]) {
                if (read_i > write_i) {
                    store[write_i] = store[read_i];
                    exists[write_i] = true;
                    exists[read_i] = false;
                }
                write_i++;
            }
            read_i++;
        }
        long amt_elements = write_i - from; // total # of elements between from and to
        return amt_elements;
    }

    //fr inclusive, to exclusive
    void spread(long from, long to, long level/*for smart rebalancing TODO*/, long amt_elements){
        long nodes = (to - from) / node_size;
        long node_fill = amt_elements / nodes;
        long plus1_nodes = amt_elements % nodes;
        long read_ptr = from + amt_elements - 1;
        for (long i = 0; i < nodes; i++){
            long my_node_fill = (i < plus1_nodes) ? node_fill+1 : node_fill;
            for (long j = 0; j < my_node_fill; j++){
                long write_ptr = from + (node_size*(nodes-i)) - (j + 1);
                store[write_ptr] = store[read_ptr];
                exists[write_ptr] = true;
                exists[read_ptr] = false;
                read_ptr--;
            }
        }
    }

    //fr inclusive, to exclusive
    void spread_spline(long from, long to, long level/*for smart rebalancing TODO*/, long amt_elements){

        long nodes = (to - from) / node_size;
        long node_fill = amt_elements / nodes;
        long plus1_nodes = amt_elements % nodes;
        long read_ptr = from + amt_elements - 1;

        auto spl = spline.upper_bound(store[read_ptr].key);
        assert(spl != spline.begin());
        spl--;
        
        auto spl_seg = spline_segments.find(spl->second.fromY);
        assert(spl_seg->first == spl->second.fromY);
        std::vector<std::pair<long, KeyType>> spl_seg_reinsert;

        for (long i = 0; i < nodes; i++){
            long my_node_fill = (i < plus1_nodes) ? node_fill+1 : node_fill;
            for (long j = 0; j < my_node_fill; j++){
                long write_ptr = from + (node_size*(nodes-i)) - (j + 1);

                if(store[read_ptr].key == spl->first  && !(write_ptr == from) && spl != spline.begin()){
                    assert(spl->second.fromY == spl_seg->first);

                    auto tmp_spl_seg = spl_seg;
                    spl_seg = std::prev(spl_seg);
                    spline_segments.erase(tmp_spl_seg);
                    spl_seg_reinsert.push_back(std::pair(write_ptr, tmp_spl_seg->second));
                        
                    spl->second.fromY = write_ptr;
                    spl--;

                    if(spl->second.toY < to){
                        auto nx = std::next(spl);
                        spl->second.toY = nx->second.fromY-1;
                        assert(nx != spline.end());
                    }                    
                }
                store[write_ptr] = store[read_ptr];
                exists[write_ptr] = true;
                exists[read_ptr] = false;
                read_ptr--;
            }
        }
        
        for (auto x : spl_seg_reinsert){
            spline_segments.emplace(x);
        }
    }

    void expand(){
        long old_sz = store.size();
        exists.resize(2*old_sz);
        store.resize(2*old_sz);
        long amt = pack(0,old_sz);

        spread(0, store.size(), 0/*unused*/, amt);

        levels++;
    }
    
    void rebalance(long pos){
        long node_ind;
        long fr;
        long to;
        long level = find_balanced_level(pos, &node_ind, &fr, &to);
        if (store.size() * 0.5 < container_size ){
            expand();
            updateSpline(0, store.size()-1);
        }
        if (PMA){
            if (level == 0){
                //Already balanced;
                //checkDistance(pos,pos);
                return;
            }
            else if (level == -1){
                //expand
                //DBGP(expand)
                //DBGPA(store.size())
                expand();
                updateSpline(0, store.size()-1);
            }
            else{
                //rebalancing
                long amt = pack(fr,to);
                DBGPB(fr, to)
                spread_spline(fr, to, level, amt);
                checkDistance(fr,to-1);
            }
        }
    }

    //Recalculates Spline in a segment
    void updateSpline(long from, long to){
        DBGP(updateSpline)
        DBGPB(from, to)

        typename std::map<long, KeyType>::const_iterator itFrom_seg = spline_segments.upper_bound(from);
        if(itFrom_seg != spline_segments.begin()){
            itFrom_seg--;
        }
        else{
            assert(false);
        }
        long fromLoc = itFrom_seg->first;
        typename std::map<KeyType, SplineSegment>::const_iterator itFrom = spline.find(itFrom_seg->second);
        assert(itFrom->second.fromY == fromLoc);

        typename std::map<long, KeyType>::const_iterator itTo_seg = spline_segments.upper_bound(to);
        if(itTo_seg != spline_segments.begin()){
            itTo_seg--;
        }
        typename std::map<KeyType, SplineSegment>::const_iterator itTo = spline.find(itTo_seg->second);
        long toLoc = itTo->second.toY;

        if (toLoc < to){
            toLoc = to;
        }
        assert(to < store.size());
            
        //DBGP(updating spline)
        //DBGPB(from, to)
        //DBGPB(itFrom->second.fromY, itFrom->second.toY)
        //DBGPB(itTo->second.fromY, itTo->second.toY)

        //erase all elements in [itFrom,itTo]
        spline_segments.erase(itFrom_seg, std::next(itTo_seg));
        spline.erase(itFrom, std::next(itTo));

        int obtained_segments = 0;
        while(fromLoc <= toLoc){
            obtained_segments++;
            SplineSegment seg;
            KeyType segMin;
            long oldFromLoc = fromLoc;
            fromLoc = getSegment(fromLoc, toLoc, EPSILON, &seg, &segMin);

            auto retPair = spline.emplace(segMin,seg);
            assert(retPair.second); //Insertion took place
            auto retPair2 = spline_segments.emplace(oldFromLoc, segMin);
            assert(retPair2.second); //Insertion took place
        }
        //std::cout << "obtained_segments= "<< obtained_segments<<std::endl;
        if (obtained_segments <= 1){
            //std::raise(SIGINT);
        }
    }
    
    inline long findNextFilled(long curr){
        while(curr < (long)exists.size() && exists[curr] == false)
        {
            curr++;
        }
        if (curr == (long)exists.size())
        {
            curr++;
        }
        return curr;
    }    

    inline double calculateDistance(double e1_x, double e1_y, double e2_x, double e2_y, double p_x, double p_y, double p1_x, double p1_y, double * sl, double * off){
        double dx = e2_x - e1_x;
        double dy = e2_y - e1_y;
        *sl = dy/dx;
        *off = e1_y - (e1_x * (*sl));

        double test_pos = p_y;
        double test_eval = (p_x * (*sl)) + (*off);
        double e = test_pos - test_eval;
        
        double test_pos2 = p1_y;
        double test_eval2 = (p1_x * (*sl)) + (*off);
        double e2 = test_pos2 - test_eval2;

        //vizSlopes(*sl, *off, e, e2);
        return std::abs(e) > std::abs(e2) ? e : e2;
    }
    
    void getSlope(double e1_x, double e1_y, double e2_x, double e2_y, double * sl, double * off){
        double dx = e2_x - e1_x;
        double dy = e2_y - e1_y;
        *sl = dy/dx;
        *off = e1_y - (e1_x * (*sl));
    }

    double getError(double p_x, double p_y, double sl, double off){
        double test_pos = p_y;
        double test_eval = (p_x * sl) + off;
        return test_pos - test_eval;
    }

    long getSegment(long start, long end, float max_epsilon, SplineSegment * seg, KeyType * segMin){
        DBGP(getSegment)
        DBGPB(start, end+1)

        std::vector<long> upper_hull;
        std::vector<long> lower_hull;
        
        long upper_anti = 0;
        long upper_anti_plus = 0;
        long lower_anti = 0;
        long lower_anti_plus = 0;
        long upper_is_edge = true;  //true -> upper is the antipodal edge
                                    //false -> lower is the antipodal edge
        end++;
        seg->fromY = start;
        while (!exists[start] && start < end){
            start++;
            if(start >= end){
                assert(false);
                seg->toY = end-1;
                return end;
            }
        }
        long i = start;

        *segMin = store[start].key;

        double slope = 0;
        double offset = start;
        double epsilon = INFINITY;
        double prevEpsilon = -INFINITY;
        double absolute_max_error = INFINITY;
                
        while(i < end)
        {
            while (!exists[i] && i < end){
                i++;
            }
            if(i >= end){
                break;
            }
            //Handle convex hull
            while (upper_hull.size() > 1){
                double cross_upper = cross((double)store[upper_hull[upper_hull.size()-2]].key, (double)upper_hull[upper_hull.size()-2], (double)store[upper_hull[upper_hull.size()-1]].key, (double)upper_hull[upper_hull.size()-1], (double)store[i].key, (double)i);
                if(cross_upper >= (double)0)
                {
                    upper_hull.pop_back();
                }
                else{
                    break;
                }
            }
            upper_hull.push_back(i);
            
            //Handle convex hull
            while (lower_hull.size() > 1){
                double cross_lower = cross((double)store[lower_hull[lower_hull.size()-2]].key, (double)lower_hull[lower_hull.size()-2], (double)store[lower_hull[lower_hull.size()-1]].key, (double)lower_hull[lower_hull.size()-1], (double)store[i].key, (double)i);
                if(cross_lower <= (double)0){
                    lower_hull.pop_back();
                }
                else{
                    break;
                }
            }
            lower_hull.push_back(i);

            //Change antipodal point position
            if (upper_hull.size() > 1 && lower_hull.size() > 1){

                if (upper_anti > upper_hull.size() - 2)
                {
                    upper_anti = upper_hull.size() - 2;
                }
                if (lower_anti > lower_hull.size() - 2)
                {
                    lower_anti = lower_hull.size() - 2;
                }
                if (upper_anti_plus > upper_hull.size() - 2)
                {
                    upper_anti_plus = upper_hull.size() - 2;
                }
                if (lower_anti_plus > lower_hull.size() - 2)
                {
                    lower_anti_plus = lower_hull.size() - 2;
                }

                if (upper_hull.size() == 2){
                    upper_anti = 0;
                    upper_anti_plus = 1;
                }
                if (lower_hull.size() == 2){
                    lower_anti = 0;
                    lower_anti_plus = 1;
                }

                {
                    if (epsilon != INFINITY){
                        prevEpsilon = epsilon;
                    }

                    {
                        double l;
                        double l_p;
                        if(lower_anti_plus == lower_hull.size()-1){
                            l = store[lower_hull[lower_anti_plus]].key - store[lower_hull[lower_anti_plus-1]].key;
                            l_p = lower_hull[lower_anti_plus] - lower_hull[lower_anti_plus-1];
                        }
                        else{
                            l = store[lower_hull[lower_anti_plus+1]].key - store[lower_hull[lower_anti_plus]].key;
                            l_p = lower_hull[lower_anti_plus+1] - lower_hull[lower_anti_plus];
                        }
                        double u = store[upper_hull[upper_anti_plus]].key - store[upper_hull[upper_anti_plus-1]].key;
                        double u_p = upper_hull[upper_anti_plus] - upper_hull[upper_anti_plus-1];
                        while (cross(u,u_p,l,l_p) < 0 && lower_anti_plus < lower_hull.size()-2){
                            lower_anti_plus++;
                            l = store[lower_hull[lower_anti_plus+1]].key - store[lower_hull[lower_anti_plus]].key;
                            l_p = lower_hull[lower_anti_plus+1] - lower_hull[lower_anti_plus];
                        }
                    }
                    {
                        double l = store[lower_hull[lower_anti_plus]].key - store[lower_hull[lower_anti_plus-1]].key;
                        double l_p = lower_hull[lower_anti_plus] - lower_hull[lower_anti_plus-1];
                        double u;
                        double u_p;
                        if(upper_anti_plus == upper_hull.size()-1){
                            u = store[upper_hull[upper_anti_plus]].key - store[upper_hull[upper_anti_plus-1]].key;
                            u_p = upper_hull[upper_anti_plus] - upper_hull[upper_anti_plus-1];
                        }
                        else{
                            u = store[upper_hull[upper_anti_plus+1]].key - store[upper_hull[upper_anti_plus]].key;
                            u_p = upper_hull[upper_anti_plus+1] - upper_hull[upper_anti_plus];
                        }
                        
                        while (cross(u,u_p,l,l_p) < 0 && upper_anti_plus < upper_hull.size()-2){
                            upper_anti_plus++;
                            u = store[upper_hull[upper_anti_plus+1]].key - store[upper_hull[upper_anti_plus]].key;
                            u_p = upper_hull[upper_anti_plus+1] - upper_hull[upper_anti_plus];
                        }                            
                    }
                    
 
                    {
                        double l = store[lower_hull[lower_anti+1]].key - store[lower_hull[lower_anti]].key;
                        double l_p = lower_hull[lower_anti+1] - lower_hull[lower_anti];
                        double u;
                        double u_p;
                        if(upper_anti_plus == upper_hull.size()-1){
                            u = store[upper_hull[upper_anti_plus]].key - store[upper_hull[upper_anti_plus-1]].key;
                            u_p = upper_hull[upper_anti_plus] - upper_hull[upper_anti_plus-1];
                        }
                        else{
                            u = store[upper_hull[upper_anti_plus+1]].key - store[upper_hull[upper_anti_plus]].key;
                            u_p = upper_hull[upper_anti_plus+1] - upper_hull[upper_anti_plus];
                        }
                        while (cross(u,u_p,l,l_p) < 0 && lower_anti < lower_hull.size()-2){
                            lower_anti++;
                            l = store[lower_hull[lower_anti+1]].key - store[lower_hull[lower_anti]].key;
                            l_p = lower_hull[lower_anti+1] - lower_hull[lower_anti];
                        }
                    }
                    {
                        double l;
                        double l_p;
                        if(lower_anti_plus == lower_hull.size()-1){
                            l = store[lower_hull[lower_anti_plus]].key - store[lower_hull[lower_anti_plus-1]].key;
                            l_p = lower_hull[lower_anti_plus] - lower_hull[lower_anti_plus-1];
                        }
                        else{
                            l = store[lower_hull[lower_anti_plus+1]].key - store[lower_hull[lower_anti_plus]].key;
                            l_p = lower_hull[lower_anti_plus+1] - lower_hull[lower_anti_plus];
                        }

                        double u = store[upper_hull[upper_anti+1]].key - store[upper_hull[upper_anti]].key;
                        double u_p = upper_hull[upper_anti+1] - upper_hull[upper_anti];
                        while (cross(u,u_p,l,l_p) < 0 && upper_anti < upper_hull.size()-2){
                            upper_anti++;
                            u = store[upper_hull[upper_anti+1]].key - store[upper_hull[upper_anti]].key;
                            u_p = upper_hull[upper_anti+1] - upper_hull[upper_anti];
                        }                            
                    }

                    //custom comparator that only compares abs value
                    auto abs_compare_lt = [&](const double& a, const double& b){
                        return std::abs(a) < std::abs(b);
                    };

                    double sl1, sl2, sl3, sl4;
                    double off1, off2, off3, off4;

                    //Diese conditions funktionieren sind aber nicht bewiesen
                    double e1 = INFINITY;
                    if (upper_anti < upper_hull.size()-1){
                        getSlope(store[upper_hull[upper_anti]].key, upper_hull[upper_anti] , store[upper_hull[upper_anti+1]].key, upper_hull[upper_anti+1], &sl1, &off1);
                        double a = getError(store[lower_hull[lower_anti]].key, lower_hull[lower_anti], sl1, off1);
                        double b = getError(store[lower_hull[lower_anti_plus]].key, lower_hull[lower_anti_plus], sl1, off1);
                        double c = getError(store[lower_hull[0]].key, lower_hull[0], sl1, off1);
                        double d = getError(store[lower_hull[lower_hull.size()-1]].key, lower_hull[lower_hull.size()-1], sl1, off1);
                        e1 = std::max({a,b,c,d}, abs_compare_lt);
                    }

                    double e2 = INFINITY;
                    if (lower_anti < lower_hull.size()-1){
                        getSlope(store[lower_hull[lower_anti]].key, lower_hull[lower_anti] , store[lower_hull[lower_anti+1]].key, lower_hull[lower_anti+1], &sl2, &off2);
                        double a = getError(store[upper_hull[upper_anti]].key, upper_hull[upper_anti], sl2, off2);
                        double b = getError(store[upper_hull[upper_anti_plus]].key, upper_hull[upper_anti_plus], sl2, off2);
                        double c = getError(store[upper_hull[0]].key, upper_hull[0], sl2, off2);
                        double d = getError(store[upper_hull[upper_hull.size()-1]].key, upper_hull[upper_hull.size()-1], sl2, off2);
                        e2 = std::max({a,b,c,d}, abs_compare_lt);
                    }

                    double e3 = INFINITY;
                    if (upper_anti > 0){
                        getSlope(store[upper_hull[upper_anti_plus-1]].key, upper_hull[upper_anti_plus-1] , store[upper_hull[upper_anti_plus]].key, upper_hull[upper_anti_plus], &sl3, &off3);
                        double a = getError(store[lower_hull[lower_anti]].key, lower_hull[lower_anti], sl3, off3);
                        double b = getError(store[lower_hull[lower_anti_plus]].key, lower_hull[lower_anti_plus], sl3, off3);
                        double c = getError(store[lower_hull[0]].key, lower_hull[0], sl3, off3);
                        double d = getError(store[lower_hull[lower_hull.size()-1]].key, lower_hull[lower_hull.size()-1], sl3, off3);
                        e3 = std::max({a,b,c,d}, abs_compare_lt);
                    }
                    double e4 = INFINITY;
                    if (lower_anti > 0){
                        getSlope(store[lower_hull[lower_anti_plus-1]].key, lower_hull[lower_anti_plus-1] , store[lower_hull[lower_anti_plus]].key, lower_hull[lower_anti_plus], &sl4, &off4);
                        double a = getError(store[upper_hull[upper_anti]].key, upper_hull[upper_anti], sl4, off4);
                        double b = getError(store[upper_hull[upper_anti_plus]].key, upper_hull[upper_anti_plus], sl4, off4);
                        double c = getError(store[upper_hull[0]].key, upper_hull[0], sl4, off4);
                        double d = getError(store[upper_hull[upper_hull.size()-1]].key, upper_hull[upper_hull.size()-1], sl4, off4);
                        e4 = std::max({a,b,c,d}, abs_compare_lt);
                    }


                    double e1_max = std::abs(e1);
                    double e2_max = std::abs(e2);
                    double e3_max = std::abs(e3);
                    double e4_max = std::abs(e4);
                
                    if (e1_max < e2_max && e1_max < e3_max && e1_max < e4_max){
                        upper_is_edge = true;
                        slope = sl1;
                        offset = off1;
                        epsilon = e1_max;
                        absolute_max_error = e1;
                    }
                    else if(e2_max < e3_max && e2_max < e4_max){
                        upper_is_edge = false;
                        slope = sl2;
                        offset = off2;
                        epsilon = e2_max;
                        absolute_max_error = e2;
                    }
                    else if(e3_max < e4_max){
                        upper_is_edge = true;
                        slope = sl3;
                        offset = off3;
                        epsilon = e3_max;
                        absolute_max_error = e3;
                    }
                    else{
                        upper_is_edge = false;
                        slope = sl4;
                        offset = off4;
                        epsilon = e4_max;
                        absolute_max_error = e4;
                    }
                        
                    if(epsilon+0.001 <= prevEpsilon){
                        DBGP(ERRRRRROR)
                        DBGP(Upper and lower hull)
                        for (auto x : upper_hull){
                            std::cout << store[x].key << " " << x << " U\n";
                        }
                        for (auto x : lower_hull){
                            std::cout << store[x].key << " " << x << " L\n";
                        }
                        vizSlopes(sl1, off1, e1, e1, store[lower_hull[lower_hull.size()-1]].key);
                        vizSlopes(sl2, off2, e2, e2, store[lower_hull[lower_hull.size()-1]].key);
                        vizSlopes(sl3, off3, e3, e3, store[lower_hull[lower_hull.size()-1]].key);
                        vizSlopes(sl4, off4, e4, e4, store[lower_hull[lower_hull.size()-1]].key);
                        std::cout << "upper_anti " << upper_anti << " lower_anti " << lower_anti << std::endl;
                        std::cout << "upper_anti_plus " << upper_anti_plus << " lower_anti_plus " << lower_anti_plus << std::endl;
                        std::cout << "epsilon " << epsilon << " prevEpsilon " << prevEpsilon << std::endl;
                        //Print Convex Hull
                        assert(false);
                    }
                    if (epsilon > max_epsilon){
                        seg->toY = i-1;
                        return i;
                    }
                }
            }
            else{
                    slope = 0;
                    offset = upper_hull[0];
                    epsilon = 0;
                    absolute_max_error = 0;
            }                
            seg->a = slope;
            seg->b = offset + (absolute_max_error/2);
            seg->e = std::abs(absolute_max_error/2);
            prevEpsilon = epsilon;

            i++;
        }
        seg->toY = end-1;
        return end;
    }


    //Debug print functions

    void assert_print(bool assertion){
        if (!assertion){
            print();
            assert(false);
        }
    }

    bool intact(){
        #ifdef NDEBUG

        return true;

        #else
        if(spline_segments.size() != spline.size()){
            assert(false);
            return false;
        }

        KeyType prev;
        bool found = false;
        for (long i = 0; i < (long)store.size(); i++){
            if (exists[i])
            {
                if (found){
                    if(store[i].key <= prev){
                        DBGP(Out of order)
                        DBGPC(store[i].key , i , prev)
                        return false;
                    }
                }
                prev = store[i].key;
                found = true;
            }
        }
        long last_to = -1;
        for (auto const& x : spline)
        {
            if (store[findNextFilled(x.second.fromY)].key != x.first){
                DBGP(bad spline indexes)
                DBGPC(store[findNextFilled(x.second.fromY)].key, x.first, x.second.fromY)
                return false;
            }
            assert(x.second.fromY<=x.second.toY);
            if (x.second.fromY > x.second.toY){
                return false;
            }
            if (last_to+1 != x.second.fromY){
                DBGP(Bad Spline)
                return false;
            }
            last_to = x.second.toY;
        }
        return true;
        #endif
    }

    bool intact_distance(){
        int num_elements = 0;
        for (int i = 0; i < store.size(); i++){
            if (exists[i]){
                KeyType val = store[i].key;
                double pos = i;

                typename std::map<KeyType, SplineSegment>::const_iterator it = spline.upper_bound(val);
                if (it != spline.begin()){
                    it--;
                }
                SplineSegment spline_seg = it->second;
                double loc = spline_seg.a * (double)val + spline_seg.b;
                double dist = pos - loc;
                if ((long)std::abs(dist) > MAX_ERROR){
                    DBGPB(pos, loc)
                    DBGPB(spline_seg.a, spline_seg.b)
                    return false;
                }
                num_elements++;
            }
        }
        if(num_elements != container_size){
            DBGPB(num_elements, container_size)
            return false;
        }
        return true;
    }

    void print_vis(){
        std::ofstream outfileData;
        outfileData.open("./viz/tmp/data.csv");

        std::cout << "X Y\n";
        outfileData << "X Y\n";
        for (long i = 0; i < (long)store.size(); i++){
            if (exists[i])
            {
                std::cout << store[i].key << " " << i << "\n";
                outfileData << store[i].key << " " << i << "\n";
            }
        }
        std::cout << "\n";
        outfileData << "\n";
        outfileData.close();
        
        std::ofstream outfileSegments;
        outfileSegments.open("./viz/tmp/segments.csv");

        std::cout << "X Y X1 Y1\n";
        outfileSegments << "X Y X1 Y1\n";

        for (auto const& x : spline)
        {
            std::cout << ((double)x.second.fromY-x.second.b)/x.second.a << " " << x.second.fromY << " " << ((double)x.second.toY-x.second.b)/x.second.a << " "  << x.second.toY <<  "\n";

            outfileSegments << ((double)x.second.fromY-x.second.b)/x.second.a << " " << x.second.fromY << " " << ((double)x.second.toY-x.second.b)/x.second.a << " "  << x.second.toY <<  "\n";
        }
        std::cout << "\n";
        outfileSegments << "\n";
        outfileSegments.close();


        std::ofstream outfileHulls;
        outfileHulls.open("./viz/tmp/hulls.csv");

        std::cout << "X Y C\n";
        outfileHulls << "X Y C\n";

        for (auto const& x : spline)
        {
            for (auto const& uh : x.second.upperHull){
                std::cout << uh.second.x << " " << uh.first << " " << "U" <<  "\n";
                outfileHulls << uh.second.x << " " << uh.first << " " << "U" <<  "\n";
            }
            for (auto const& lh : x.second.lowerHull){
                std::cout << lh.second.x << " " << lh.first << " " << "L" <<  "\n";
                outfileHulls << lh.second.x << " " << lh.first << " " << "L" <<  "\n";
            }
        }
        std::cout << "\n";
        outfileHulls << "\n";
        outfileHulls.close();
    }

    void print(){
        DBGPA(levels);
        for (long i = 0; i < (long)store.size(); i++){
            if (exists[i])
            {
                std::cout << i << "->" << store[i].key << " ";
            }
            else
            {
                std::cout << "- ";
            }
        }
        std::cout << "\n";

        std::cout << "Spline \n";
        for (auto const& x : spline)
        {
            std::cout << x.first << "->" << x.second.fromY << ":" << x.second.toY << std::endl;
        }
        std::cout << "\n";
    }

    void vizSlopes(double sl, double off, double err, double err2){
        double le = 1000;
        std::cout << 0 << " " << off << " " << le << " " << (le*sl) + off << "\n";
        std::cout << 0 << " " << off + err << " " << le << " " << (le*sl) + off + err << "\n";
        std::cout << 0 << " " << off + err2 << " " << le << " " << (le*sl) + off + err2 << "\n";
    }

    void debugSegments(){
        debugSegments(0,store.size(), true);
    }

    void debugSegments(long fr, long to, bool sparse){
        double sum_dists = 0;
        long num_elements = 0;
        for (int i = fr; i < to; i++){
            if (exists[i]){
                KeyType val = store[i].key;
                double pos = i;

                typename std::map<KeyType, SplineSegment>::const_iterator it = spline.upper_bound(val);
                if (it != spline.begin()){
                    it--;
                }
                SplineSegment spline_seg = it->second;
                double loc = spline_seg.a * (double)val + spline_seg.b;
                double dist = pos - loc;
                if ((long)std::abs(dist) > MAX_ERROR){
                    std::cout << "\033[0;31m" << val << " : " << dist << " loc = " << (long)loc << "\033[0m \n";
                }
                else{
                    if(!sparse){
                        std::cout << val << " : " << dist << "\n";
                    }
                }
                num_elements++;
                sum_dists += std::abs(dist);
            }
        }
        std::cout << "average distance = " << sum_dists/num_elements << "\n";
        std::cout << "num elements = " << num_elements << "\n";
        std::cout << "wanted container size = " << container_size << "\n";
        std::cout << "num segments = " << spline.size() << "\n";
        std::cout << "average segment size = " << num_elements / spline.size() << "\n";
    }

    private:
};