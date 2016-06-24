// Timing_Stats.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// store and report timing performance statistics for any routine

#ifndef TIMING_STATS_HPP
#define TIMING_STATS_HPP

#include <time.h>

class Timing_Stats
{
  clock_t start_time, end_time;
public:
  double cpu_time;
  Timing_Stats() : start_time(clock()), cpu_time(0.) {}
  virtual void start_clock()
  {
    start_time = clock();
  }
  virtual double collect_stats()
  {
    end_time = clock();
    const double time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    cpu_time += time;
    start_time = clock();
    return time;
  }
};

class Wheel_Stats : public Timing_Stats
{
public:
  double cpu_time_success;
  size_t spins, spins_success;
  Wheel_Stats() : Timing_Stats(), cpu_time_success(0), spins(0), spins_success(0) {}
  virtual double collect_stats(bool valid_dart)
  {
    double time = Timing_Stats::collect_stats();
    ++spins;
    if (valid_dart)
    {
      ++spins_success;
      cpu_time_success += time;
    }
    return time;
  }
};

class Rebalance_Stats : public Timing_Stats
{
public:
  Rebalance_Stats() : Timing_Stats(), num_rebalances(0) {}
  size_t num_rebalances;
};

#endif
