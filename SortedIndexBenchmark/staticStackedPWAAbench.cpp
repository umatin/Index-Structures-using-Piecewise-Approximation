#define NDEBUG
#include "benchmark.h"
#include "../StackedPwApproximation/stackedSpline.h"
#include <string>

int main(int argc, char **argv) {
    if(argc < 3){
        return 0;
    }
    std::string path(argv[1]);

    Benchmark<u_int64_t> b;
    stackedSpline<u_int64_t,u_int64_t> st(strtol(argv[2], nullptr, 0));

    b.loadData(path);
    b.sortData();
    b.sampleLookups(200000000);

    u_int64_t timing_insert = b.Build(st);
    u_int64_t timing_lookup = b.Lookup(st);

    std::cout << "timing_insert = " << b.nano_seconds_to_double(timing_insert) << std::endl;
    std::cout << "timing_lookup = " << b.nano_seconds_to_double(timing_lookup) << std::endl;

    st.print();
}