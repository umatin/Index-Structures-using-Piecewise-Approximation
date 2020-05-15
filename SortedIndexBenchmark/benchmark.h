
#include <vector>
#include <algorithm>
#include <random>

#include "util.h"

template<typename KeyType = uint64_t>
class Benchmark{
    std::vector<KeyType> data_;
    std::vector<KeyType> lookups_;
public:
    Benchmark(){

    }

    void loadData(const std::string& filename){
        std::cout << "Loading data: " << filename << std::endl;
        data_ = load_data<KeyType>(filename);
    }
    
    void sortData(){
        std::sort(data_.begin(), data_.end());
    }
    
    void shuffleData(){
        std::random_shuffle(data_.begin(), data_.end());
    }

    void sampleData(long num){
        std::cout << "Sampling " << num << " from data" << std::endl;
        assert(data_.size() > 0);
        
        std::vector<KeyType> data_2;

        for(int i = 0; i < num; i++){
            long loc = rand() % data_.size();
            data_2.push_back(data_[loc]);
        }
        
        data_ = data_2;
    }

    void makeRandomData(long num){
        std::cout << "Loading data: Random Data" << num << std::endl;
        data_.clear();
        for(int i = 0; i < num; i++){
            data_.push_back(rand());
        }
    }

    void loadLookups(const std::string& filename){
        std::cout << "Loading lookups: " << filename << std::endl;
        lookups_ = load_data<KeyType>(filename);
    }

    void sampleLookups(long num){
        std::cout << "Sampling " << num << " from data" << std::endl;
        assert(data_.size() > 0);

        lookups_.clear();
        for(int i = 0; i < num; i++){
            long loc = rand() % data_.size();
            lookups_.push_back(data_[loc]);
        }
    }
    template<class Index>
    uint64_t Build(Index & index_) {
        auto start = std::chrono::high_resolution_clock::now();
        index_.build(data_);
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    template<class Index, bool debug = false>
    uint64_t Insert(Index & index_) {
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < data_.size(); i++){
            KeyType key = data_[i];
            index_.insert(key,key);
            
            if(debug){
                if(i % 1000000 == 0 && i != 0){
                    std::cout << "i->" << i << std::endl;
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    template<class Index>
    std::vector<uint64_t> Insert_timeEach(Index & index_) {
        std::vector<uint64_t> timings;
        for (int i = 0; i < data_.size(); i++){
            KeyType key = data_[i];
            auto start = std::chrono::high_resolution_clock::now();
            index_.insert(key,key); 
            auto end = std::chrono::high_resolution_clock::now();
            timings.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        }
        return timings;
    }

    template<class Index>
    uint64_t Lookup(Index & index_) {
        std::cout << "timing " << lookups_.size() << " lookups" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < lookups_.size(); i++){
            KeyType key = lookups_[i];
            KeyType val = index_.get(key);
            if(key != val){
                std::cout << "key=" << key << " != val=" << val << std::endl;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    template<class Index>
    std::vector<uint64_t> Lookup_timeEach(Index & index_) {
        std::cout << "timing " << lookups_.size() << " lookups each" << std::endl;
        std::vector<uint64_t> timings;
        for (int i = 0; i < lookups_.size(); i++){
            
            KeyType key = lookups_[i];
            auto start = std::chrono::high_resolution_clock::now();
            KeyType val = index_.get(key);
            auto end = std::chrono::high_resolution_clock::now();
            timings.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

            if(key != val){
                std::cout << "key=" << key << " != val=" << val << std::endl;
            }
        }
        return timings;
    }

    static double nano_seconds_to_double(u_int64_t t){
        std::chrono::nanoseconds t_ns(t);

        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_ns);
        return time_span.count();
    }

};