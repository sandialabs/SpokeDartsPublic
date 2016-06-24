//  Search_Stats.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Search_Stats.hpp"

void Search_Stats::report_stats(std::ostream &out)
{
  out << "Search Stats: " << _searches << " searches in " << cpu_time << " seconds." << std::endl;
  out << "  found " << _found << " of " << _visited << " visited of " << _searched << " searched." << std::endl;

  const double num_neighbors = _found / ((double) _searches );
  out << "  Average number found per search = " << num_neighbors << std::endl;

  const double tff = _found / ( (double) _searched );
  const double tvf = _visited / ( (double) _searched );
  const double tvff = _found / ( (double) _visited );
  out << "  Total fraction visited/total = " << tvf << ", found/total = " << tff << ", found/visited = " << tvff << std::endl;

  const double ff = _found_fraction / _searches;
  const double vf = _visited_fraction / _searches;
  const double vff = _visited_found_fraction / _searches;
  out << "  Average fraction visited/total = " << vf << ", found/total = " << ff << ", found/visited = " << vff << std::endl;
  
}
