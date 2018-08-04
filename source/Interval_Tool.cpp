// Interval_Tool.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Interval_Tool.hpp"
#include <assert.h>
#include <cstddef>

void Interval_Tool::pick_interval_uniformly( double t, const std::vector<double> &intervals, double sum_lengths, double &A, double &B )
{
  const size_t num_intervals = intervals.size() / 2;
  double seg_sum = 0;
  assert( t <= 1. && t > 0);
  double seg_sum_limit = t * sum_lengths;
  for (size_t i = 0; i < num_intervals; ++i)
  {
    A = intervals[i*2];
    B = intervals[i*2+1];
    assert( B >= A );
    seg_sum += B-A;
    if ( seg_sum > seg_sum_limit )
      return;
  }
}

void Interval_Tool::longest_interval( const std::vector<double> &intervals, double &A, double &B )
{
  const size_t num_intervals = intervals.size() / 2;
  A = 0.; B = 0.;
  for (size_t i = 0; i < num_intervals; ++i)
  {
    double a = intervals[i*2];
    double b = intervals[i*2+1];
    assert( b >= a );
    if (b-a > B-A)
    {
      B = b;
      A = a;
    }
  }
}
