#include <chrono>
#include <iostream>
#include <ratio>
#include <ctime>

int main()
{
    using namespace std::chrono;
    typedef std::chrono::seconds sec;

    steady_clock::time_point start = steady_clock::now();

    std::cout << "printing 1000 ints need: \n";
    for (int i = 0; i < 1000; ++i) {
        std::cout << i << " ";
    }
    std::cout << "\n";
    steady_clock::time_point stop = steady_clock::now();
    
    duration<double> time_span = duration_cast<duration<double>>(stop-start);
    std::cout << "It tooks: " << time_span.count() << "seconds\n";
    
    return 0;

}
