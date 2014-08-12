#ifndef SIMULATORTIMER_HEADER_INCLUDED
#define SIMULATORTIMER_HEADER_INCLUDED

#include <iosfwd>
#include <vector>

class SimulatorTimer
{
public:
    /// Default constructor.
    SimulatorTimer();

    /// Initialize from parameters. Accepts the following:
    ///    num_psteps    (default 1)
    ///    stepsize_days (default 1)
    void init(const double nsteps, const double dt);

    /// Total number of steps.
    int numSteps() const;
    /// Current step number. This is the number of timesteps that
    /// has been completed from the start of the run. The time
    /// after initialization but before the simulation has started
    /// is timestep number zero.
    int currentStepNum() const;
    /// Set current step number.
    void setCurrentStepNum(int step);

    /// Current step length. This is the length of the step
    /// the simulator will take in the next iteration.
    /// @note if done(), it is an error to call currentStepLength().
    double currentStepLength() const;

    /// Previous step length. This is the length of the step that
    /// was taken to arrive at this time.
    ///
    /// @note if no increments have been done (i.e. the timer is
    /// still in its constructed state and currentStepNum() == 0),
    /// it is an error to call stepLengthTaken().
    double stepLengthTaken () const;


    /// Time elapsed since the start of the simulation until the
    /// beginning of the current time step [s].
    double simulationTimeElapsed() const;


    /// Total time.
    double totalTime() const;

    /// Set total time.
    /// This is primarily intended for multi-epoch schedules,
    /// where a timer for a given epoch does not have
    /// access to later timesteps.
    void setTotalTime(double time);

    /// Print a report with current and total time etc.
    /// Note: if done(), it is an error to call report().
    void report(std::ostream& os) const;

    /// Next step.
    SimulatorTimer& operator++();

    /// Return true if op++() has been called numSteps() times.
    bool done() const;

private:
    std::vector<double> timesteps_;
    int current_step_;
    double current_time_;
    double total_time_;
};

#endif // SIMULATORTIMER_HEADER_INCLUDED
