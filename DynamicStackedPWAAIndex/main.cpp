#define NDEBUG
#include "spline_tree.h"
#include <vector>

#include <functional>
#include <chrono>

int main(){
    srand(time(NULL));
    dynamicPWAAindex<long, long> sts(4,16);
    std::vector<long> data;
    int num_inserts = 10000000;
    for(int i = 0; i < num_inserts; i++){
        long val = rand();
        sts.insert(val,val);
        data.push_back(val);
        if((i) % 1000000 == 0){
            std::cout << "i->" << i << std::endl;
        }
    }

    DBGP(doing lookups)

    auto start = std::chrono::high_resolution_clock::now();
    for(long i = 0; i < data.size(); i++){
        long ret = sts.get(data[i]);
        if (ret != data[i]){
            sts.get<false>(data[i]);
            std::cout << "errror" << std::endl;
            std::cout << "ret=" << ret << "data[i]=" << data[i] << std::endl;
            //sts.print_contains(data[i]);
            break;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << "Lookup Time: " << time_span.count() << "s" << std::endl;
    



    std::cout << "sts.count() = " << sts.count() << std::endl;
    std::cout << "sts.count_base_segments() = " << sts.count_base_segments() << std::endl;
    std::cout << "sts.get_height() = " << sts.get_height() << std::endl;
    std::cout << "total_inserts = " << total_inserts << std::endl;
}