#include <opm/lattice/SimulatorTimer.hpp>
#include <ostream>
#include <numeric>
#include <cassert>

/// Default constructor.
SimulatorTimer::SimulatorTimer()
    : current_step_(0),
      current_time_(0.0)
{
}

/// Initialize from parameters. Accepts the following:
///    num_psteps    (default 1)
///    stepsize_days (default 1)
void SimulatorTimer::init(const double nsteps, const double dt)
{
    const int num_psteps = nsteps;
    const double stepsize = dt;
    timesteps_.clear();
    timesteps_.resize(num_psteps, stepsize);
    total_time_ = num_psteps*stepsize;
}

/// Total number of steps.
int SimulatorTimer::numSteps() const
{
    return timesteps_.size();
}

/// Current step number.
int SimulatorTimer::currentStepNum() const
{
    return current_step_;
}

/// Set current step number.
void SimulatorTimer::setCurrentStepNum(int step)
{
    current_step_ = step;
    current_time_ = std::accumulate(timesteps_.begin(), timesteps_.begin() + step, 0.0);
}


/// Current step length.
double SimulatorTimer::currentStepLength() const
{
    assert(!done());
    return timesteps_[current_step_];
}

double SimulatorTimer::stepLengthTaken() const
{
    assert(current_step_ > 0);
    return timesteps_[current_step_ - 1];
}

/// time elapsed since the start of the simulation [s].
double SimulatorTimer::simulationTimeElapsed() const
{
    return current_time_;
}




/// Total time.
double SimulatorTimer::totalTime() const
{
    return total_time_;
}

    /// Set total time.
    /// This is primarily intended for multi-epoch schedules,
    /// where a timer for a given epoch does not have
    /// access to later timesteps.
void SimulatorTimer::setTotalTime(double time)
{
    total_time_ = time;
}

/// Print a report with current and total time etc.
void SimulatorTimer::report(std::ostream& os) const
{
    os << "\n\n---------------    Simulation step number " << currentStepNum() << "    ---------------"
       << "\n      Current time (days)     " << simulationTimeElapsed()
       << "\n      Current stepsize (days) " << currentStepLength()
       << "\n      Total time (days)       " << totalTime()
       << "\n" << std::endl;
    }

/// Next step.
SimulatorTimer& SimulatorTimer::operator++()
{
    assert(!done());
    current_time_ += timesteps_[current_step_];
    ++current_step_;
    return *this;
}

    /// Return true if op++() has been called numSteps() times.
bool SimulatorTimer::done() const
{
    return int(timesteps_.size()) == current_step_;
}

