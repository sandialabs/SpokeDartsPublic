// Interval_Tool.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// abstract data, nothing about the domain dimension or geometry

#ifndef INTERVAL_TOOL_HPP
#define INTERVAL_TOOL_HPP

#include <vector>

// intervals are defined by start-stop A-B distance pairs, in the vector
class Interval_Tool
{
public: 
  void pick_interval_uniformly( double t, const std::vector<double> &intervals, double sum_lengths, double &A, double &B );
  void longest_interval( const std::vector<double> &intervals, double &A, double &B );
};

#endif