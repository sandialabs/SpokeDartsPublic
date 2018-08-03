//  Search_Factory.h
//  spokes


// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef __spokes__Search_Factory__
#define __spokes__Search_Factory__

#include "Search_Structure.hpp"

class Search_Factory
{
public:
  enum Search_Type { UNSPECIFIED=0, ARRAY, GRID, TREE, RANGE };
  const std::vector<const std::string> search_description = { "unspecified", "exhaustive search array, O(n)", "a uniform background grid", "a k-d tree", "a range tree (experimental!)" };
  const std::vector<const std::string> search_name = { "unspecified", "array", "grid", "k-d tree", "range tree (experimental!)" };

  static
  Search_Structure* new_search( Search_Type search_type, Spheres *spheres, bool is_global = false, double search_distance = 0.2, double xmax = 1, double xmin = 0.);
  
  static
  void delete_search( Search_Structure *search )
  { delete search; }

};

#endif /* defined(__spokes__Search_Factory__) */
