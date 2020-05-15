#include <vector>
#include <cmath>
#include <cassert>

template<class K, class V>
class key_val{
public:
    K key;
    V val;
    key_val(K k, V v){
        this->key = k;
        this->val = v;
    }
};

template<class K, class V>
long cross(long x1_o, long x1, long x2_o, long x2, key_val<K,V> * store){
    return (((long)store[x1].key - (long)store[x1_o].key) * (x2 - x2_o)) - ((x1 - x1_o) * ((long)store[x2].key - (long)store[x2_o].key));
}

template<class K, class V>
long cross(long o, long x1, long x2, key_val<K,V> * store){
    return (((long)store[x1].key - (long)store[o].key) * (x2 - o)) - ((x1 - o) * ((long)store[x2].key - (long)store[o].key));
}

class SplineSegment{
public:
    double a;
    double b;
    long predict(double x){
        return (long) ((a*x)+b);
    }
    SplineSegment(double a, double b){
        this->a = a;
        this->b = b;
    }
    SplineSegment(){
    }
};


template<class K, class V>
inline double calculateDistance(long e1, long e2, long p, key_val<K,V> * store, double * a, double * b){
        double dx = (double)store[e2].key - (double)store[e1].key;
        double dy = (double)e2 - (double)e1;
        * a = dy/dx;
        * b = e1 - ((double)store[e1].key * (* a));

        double test_pos = p;
        double test_eval = ((*a) * store[p].key) + (*b);
        double e = test_pos - test_eval;

        return e;
    }

template<class K, class V>
long getSegment(long start, long end, double max_epsilon, K * key, SplineSegment * seg, key_val<K,V> * store){
    std::vector<long> upper_hull;
    std::vector<long> lower_hull;
    
    long i = start;
    long upper_anti = 0;
    long lower_anti = 0;
    long last_e_0 = 0;
    long last_e_1 = 1;

    *key = store[start].key;

    double slope = -1;
    double offset = 0;
    double prev_epsilon = -INFINITY;
    double signed_epsilon = 0;
            
    while(i < end)
    {

        while (upper_hull.size() > 1){
            double cross_upper = cross(upper_hull[upper_hull.size()-2], upper_hull[upper_hull.size()-1], i, store);
            if(cross_upper >= 0){
                upper_hull.pop_back();
            }
            else{
                break;
            }
        }
        upper_hull.push_back(i);
        
        while (lower_hull.size() > 1){
            double cross_lower = cross(lower_hull[lower_hull.size()-2], lower_hull[lower_hull.size()-1], i, store);
            if(cross_lower <= 0){
                lower_hull.pop_back();
            }
            else{
                break;
            }
        }
        lower_hull.push_back(i);

        //Change antipodal point position
        if (upper_hull.size() > 1 && lower_hull.size() > 1){
            bool upper_changed = upper_anti > upper_hull.size() - 2;
            bool lower_changed = lower_anti > lower_hull.size() - 2;

            //assert(!upper_changed != !lower_changed);

            if (upper_changed)
            {
                upper_anti = upper_hull.size() - 2;
            }
            if (lower_changed)
            {
                lower_anti = lower_hull.size() - 2;
            }

            if(upper_changed || cross(last_e_0,last_e_1,upper_hull[upper_anti], upper_hull[upper_anti+1], store) > 0){
                last_e_0 = upper_hull[upper_anti];
                last_e_1 = upper_hull[upper_anti+1];

                while(lower_anti < lower_hull.size() - 2 && cross(last_e_0, last_e_1, lower_hull[lower_anti], lower_hull[lower_anti+1], store) < 0){
                    lower_anti++;
                }
                if(lower_hull[lower_anti] >= upper_hull[upper_anti]){
                    signed_epsilon = calculateDistance(last_e_0, last_e_1, lower_hull[lower_anti], store, &slope, &offset);
                }
                else{
                    std::cout << "Hello" << std::endl;
                    while(lower_hull[lower_anti+1] < upper_hull[upper_anti]){
                        lower_anti++;
                    }
                    last_e_0 = lower_hull[lower_anti];
                    last_e_1 = lower_hull[lower_anti+1];
                    signed_epsilon = calculateDistance(last_e_0, last_e_1, upper_hull[upper_anti], store, &slope, &offset);
                }
            }
            else if(lower_changed || cross(last_e_0,last_e_1,lower_hull[lower_anti], lower_hull[lower_anti+1], store) < 0){
                last_e_0 = lower_hull[lower_anti];
                last_e_1 = lower_hull[lower_anti+1];

                while(upper_anti < upper_hull.size() - 2 && cross(last_e_0, last_e_1, upper_hull[upper_anti], upper_hull[upper_anti+1], store) > 0){
                    upper_anti++;
                }
                if(lower_hull[lower_anti] <= upper_hull[upper_anti]){
                    signed_epsilon = calculateDistance(last_e_0, last_e_1, upper_hull[upper_anti], store, &slope, &offset);
                }
                else{
                    std::cout << "Hello" << std::endl;
                    while(upper_hull[upper_anti+1] < lower_hull[lower_anti]){
                        upper_anti++;
                    }
                    last_e_0 = upper_hull[upper_anti];
                    last_e_1 = upper_hull[upper_anti+1];
                    signed_epsilon = calculateDistance(last_e_0, last_e_1, lower_hull[lower_anti], store, &slope, &offset);
                }
            }
            /*if (std::abs(signed_epsilon) < prev_epsilon){
                std::cout << signed_epsilon << " - " << prev_epsilon << std::endl;
            }

            if (std::abs(signed_epsilon) != INFINITY){
                prev_epsilon = std::abs(signed_epsilon);
            }*/

            if (std::abs(signed_epsilon) > max_epsilon){
                return i;
            }
            
            seg->a = slope;
            seg->b = offset + (signed_epsilon/2);
        }
        i++;
    }
    if (slope == -1){
        signed_epsilon = calculateDistance(last_e_0, last_e_1, lower_hull[lower_anti], store, &slope, &offset);
        seg->a = slope;
        seg->b = offset + (signed_epsilon/2);
    }
    return i;
}

template<class K, class V>
inline long exponentialSearch(K x, int off, key_val<K,V> * d, long d_size) 
{
    if (d[off].key == x){
        return off; 
    }
    long i = 1; 
    while (off+i < d_size ? d[off+i].key <= x : false){
        i = i*2; 
    }
    return binarySearch(x, off+(i/2), std::min(off+i+1, d_size), d); 
}

template<class K, class V>
inline long exponentialSearchBW(K x, int off, key_val<K,V> * d, long d_size) 
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


template<class K, class V>
inline long binarySearch(K x, long fr, long to, key_val<K,V> * d){
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

template<class K, class V>
inline long get_pos_leq_exponentialSearch(K x, long prediction, key_val<K,V> * d, long d_size){
    long i = prediction;
    i = std::max(i,0l);
    i = std::min(i,d_size-1);
    
    if(d[i].key < x){
        long ret = exponentialSearch(x, i, d, d_size);
        ret--;
        return ret;
    }
    else if(d[i].key > x){
        long ret = exponentialSearchBW(x, i, d, d_size);
        ret--;
        return ret;
    }
    else{
        return i;
    }
}


template<class K, class V>
inline long get_pos_leq(K x, long prediction, key_val<K,V> * d, long container_size){
    long i = prediction;
    i = std::max(i,0l);
    i = std::min(i,container_size-1);
    
    long move = 0;
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
        return i;
    }
    else{
        return i;
    }
}

template<class K, class V>
inline long get_pos_leq(K x, long prediction, key_val<K,V> * d, long container_size, long * move){
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
                *move = *move+1;
            }
        }
        i--;
        return i;
    }
    else if(d[i].key > x){
        while(i > 0){  // this is intentional
            if(d[i].key <= x){
                break;
            }
            else{
                i--;
                *move = *move+1;
            }
        }
        return i;
    }
    else{
        return i;
    }
}