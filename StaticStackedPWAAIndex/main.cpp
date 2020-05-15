#include "my_util.h"
#include "stackedPWAA.h"
#include <stdlib.h>
#include <chrono>
#include <algorithm>

int main()
{
    srand(time(NULL));
    std::vector<long> d;
    std::string line;
    std::ifstream myfile ("dataset.cfg");
    int amt_data = 20000000;
    if (myfile.is_open())
    {
        getline(myfile,line);
    }

    if(line.length()>1){
        std::cout << "loading data " << line << "\n";
        myfile.close();
        d = load_data<long>(line);
        std::sort(d.begin(), d.end());
    }
    else{
        std::cout << "random Data\n";
        srand (time(NULL));
        for (int i = 0; i < amt_data; i++)
        {
            d.push_back(rand());
            //d.push_back(i);
        }
        std::sort(d.begin(), d.end());
    }
    
    stackedPWAA<long,long> s(16);

    std::cout << "Building" << std::endl;
    s.build(d);
    std::cout << "datasize: " << d.size() << std::endl;
    std::cout << "store size: " << s.store.size() << std::endl;
    std::cout << "spline size: " << s.spline_stack.size() << std::endl;
    for(int i = 0; i < s.spline_stack.size(); i++){
        std:: cout << "layer: " << i << " spl size: " << s.spline_stack[i].size() << std::endl;
    }

    std::cout << "Timing" << std::endl;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < amt_data; i++){
        long testKey = d[rand()%amt_data];
        long r = s.get(testKey);
        if (r != testKey){
            std::cout << "error\n";
        }
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Time: " << time_span.count() << "s" << std::endl;;

}