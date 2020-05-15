#define NDEBUG
#include "benchmark.h"
#include "../SplineTree/spline_tree.h"
#include <string>

int main(int argc, char **argv) {
    if(argc < 3){
        return 0;
    }
    std::string path(argv[1]);

    Benchmark<u_int64_t> b;
    spline_tree<u_int64_t,u_int64_t> st(4, strtol(argv[2], nullptr, 0));

    b.loadData(path);
    b.shuffleData();
    b.sampleLookups(200000000);

    u_int64_t timing_insert = b.Insert<spline_tree<u_int64_t,u_int64_t>, true>(st);
    u_int64_t timing_lookup = b.Lookup(st);

    std::cout << "timing_insert = " << b.nano_seconds_to_double(timing_insert) << std::endl;
    std::cout << "timing_lookup = " << b.nano_seconds_to_double(timing_lookup) << std::endl;
    
    std::cout << "st.count_base_segments " << st.count_base_segments() << std::endl;
    std::cout << "st.height " << st.height << std::endl;
}