#define NDEBUG //To disable all asserts and most debug prints


#include "my_util.h"
#include "simple_PWAA.h"
#include <stdlib.h>
#include <chrono>


KeyValue<u_int64_t> makeKV(u_int64_t v){
    return (KeyValue<u_int64_t>){v,v};
}
int main()
{
    simple_PWAA<u_int64_t, true> p;

    std::string line;
    std::ifstream myfile ("dataset.cfg");
    std::vector<u_int64_t> data;
    int amt_data = 50000000;

    if (myfile.is_open())
    {
        getline (myfile,line);
        if(line.length()>1){
            std::cout << line << '\n';
            myfile.close();
            std::vector<u_int64_t> d = load_data<u_int64_t>(line);
            for (int i = 0; i < amt_data; i++){
                auto x = d[i];
                data.push_back(x);
            }
        }
        else{
            std::cout << "random Data\n";
            //srand (time(NULL));

            for (int i = 0; i < amt_data; i++)
            {
                u_int64_t val = rand();
                data.push_back(val);
            }
        }
        std::random_shuffle(data.begin(), data.end());
        
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        for (auto val : data){
            p.insert(makeKV(val));
        }
        p.insert(makeKV(data[0]));
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Insert Time: " << time_span.count() << "s" << std::endl;

        std::random_shuffle(data.begin(), data.end());

        std::chrono::high_resolution_clock::time_point t1_lookup = std::chrono::high_resolution_clock::now();
        for (auto val : data){
            KeyValue<u_int64_t>* x = p.get(val);
            if (x->key != val){
                std::cout << "error" << val << std::endl;
            }
        }
        std::chrono::high_resolution_clock::time_point t2_lookup = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span_lookup = std::chrono::duration_cast<std::chrono::duration<double>>(t2_lookup - t1_lookup);
        std::cout << "Lookup Time: " << time_span_lookup.count() << "s" << std::endl;
    }

    
    //p.print();
    //p.updateSpline(0,p.store.size()-1);
    p.debugSegments();
    DBGP(exited normally)
}