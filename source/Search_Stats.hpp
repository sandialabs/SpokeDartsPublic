//
//  Search_Stats.hpp
//  Search_Stats.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// store and report performance statistics for proximity searches

#ifndef SPOKES_SEARCHSTATS_HPP
#define SPOKES_SEARCHSTATS_HPP

#include <iostream>
#include "Timing_Stats.hpp"

class Search_Stats : public Timing_Stats
{
public:
  Search_Stats() :
  Timing_Stats(),
  _searches(0), _searched(0), _visited(0), _found(0), _visited_fraction(0.), _found_fraction(0.), _visited_found_fraction(0.) {}
  
  virtual void collect_stats(size_t searched, size_t visited, size_t found)
  {
    ++_searches;
    _searched += searched;
    _visited += visited;
    _found += found;
    const double vf = searched ? ((double) visited) / ((double) searched) : 0.;
    _visited_fraction += vf;
    const double ff = searched ? ((double) found) / ((double) searched) : 0.;
    _found_fraction += ff;
    const double vff = visited ? ((double) found) / ((double) visited) : 0.;
    _visited_found_fraction += vff;
    
     Timing_Stats::collect_stats();
  }
  
  void report_stats(std::ostream &out = std::cout);
  
  // read only access
  size_t searches() const { return _searches; }
  size_t searched() const { return _searched; }
  size_t visited() const { return _visited; }
  size_t found() const { return _found; }
  
  const double average_neighbors() const {return _found / ((double) _searches ); }
  
protected:
  // data
  size_t _searches, _searched, _visited, _found;
  double  _visited_fraction, _found_fraction, _visited_found_fraction;
};

#endif
