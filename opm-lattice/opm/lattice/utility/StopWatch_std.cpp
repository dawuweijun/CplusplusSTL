#include <opm/lattice/StopWatch.hpp>
#include <iostream>
#include <ctime>
#include <ratio>
#include <cassert>
StopWatch::StopWatch()
    : state_(UnStarted)
{
}


void StopWatch::start()
{
    start_time_ = std::chrono::steady_clock::now();
    last_time_ = start_time_;
    state_ = Running;
}

void StopWatch::stop()
{
    if (state_ != Running) {
        std::cout << "Called stop() on a StopWatch that was not running." << std::endl;
        exit(1);
    }
    stop_time_ = std::chrono::steady_clock::now();
    state_ = Stopped;
}

double StopWatch::secsSinceLast()
{
    std::chrono::steady_clock::time_point run_time;
    if (state_ == Running) {
    	run_time = std::chrono::steady_clock::now();
    } else if (state_ == Stopped) {
    	run_time = stop_time_;
    } else {
    	assert(state_ == UnStarted);
    	std::cout << "Called secsSinceLast() on a StopWatch that had not been started.\n";
        exit(1);
    }
    
    std::chrono::duration<double> dur = std::chrono::duration_cast<std::chrono::duration<double>>(run_time - last_time_);
    last_time_ = run_time;

    return double(dur.count());
}

double StopWatch::secsSinceStart()
{
    std::chrono::steady_clock::time_point run_time;
    if (state_ == Running) {
		run_time = std::chrono::steady_clock::now();
    } else if (state_ == Stopped) {
		run_time = stop_time_;
    } else {
		assert(state_ == UnStarted);
		std::cout <<"Called secsSinceStart() on a StopWatch that had not been started.\n";
        exit(1);
    }
    std::chrono::duration<double> dur = std::chrono::duration_cast<std::chrono::duration<double>>(run_time - start_time_);
    return double(dur.count());
}
