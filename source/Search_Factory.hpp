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
  enum Search_Type { UNSPECIFIED, ARRAY, GRID, TREE, RANGE };

  static
  Search_Structure* new_search( Search_Type search_type, Spheres *spheres, bool is_global = false, double search_distance = 0.2, double xmax = 1, double xmin = 0.);
  
  static
  void delete_search( Search_Structure *search )
  { delete search; }

};

#endif /* defined(__spokes__Search_Factory__) */
