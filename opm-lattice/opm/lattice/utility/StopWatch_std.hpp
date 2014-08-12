#ifndef STOPWATCH_HEADER_INCLUDED
#define STOPWATCH_HEADER_INCLUDED

#include <chrono>

class StopWatch
{
public:
    /// Default constructor. Before the StopWatch is start()-ed,
    /// it is an error to call anything other than start().
    StopWatch();

    /// Starts the StopWatch. It is always legal to call
    /// start(), even if not stop()-ped.
    void start();
    /// Stops the StopWatch. The watch no longer runs, until
    /// restarted by a call to start().
    void stop();

    /// \return the number of running seconds that have passed
    /// since last call to start(), secsSinceLast() or
    /// secsSinceStart()
    double secsSinceLast();
    /// \return the number of running seconds that have passed
    /// since last call to start().
    double secsSinceStart();

private:
    enum StopWatchState { UnStarted, Running, Stopped };

    StopWatchState state_;
    std::chrono::steady_clock::time_point start_time_;
    std::chrono::steady_clock::time_point last_time_;
    std::chrono::steady_clock::time_point stop_time_;
};

#endif // OPM_STOPWATCH_HEADER
