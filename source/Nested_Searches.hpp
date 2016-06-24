// Nested_Searches.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef NESTED_SEARCHES
#define NESTED_SEARCHES

// A nested chain of searches
#include <iostream>
#include "Search_Structure.hpp"
#include "Sphere_Container.hpp"

class Nested_Searches
{
public:
  Nested_Searches() : _subtract_radius(false) {}

  // the first link is used to find all items within distance and the results are passed to the second link, etc.
  // The first link finds all items within the first distance
  std::vector< Search_Structure * > _searches;
  std::vector< double > _distances;
  bool _subtract_radius;

  void clear() { _searches.clear(); _distances.clear(); }
  void reserve( size_t sz ) { _searches.reserve(sz); _distances.reserve(sz-1); }

  void gather_neighbors( const double *p )
  {
    assert( _searches.size() == _distances.size() + 1);
    double old_dist = std::numeric_limits<double>::max();
    for ( size_t i = 0; i < _searches.size() - 1; ++i )
    {
      assert( _distances[i] <= old_dist );
 
      _searches[i]->all_near_spheres( _searches[i+1]->sphere_container(), p, _distances[i], _subtract_radius );

      Sphere_Container *sc = _searches[i]->sphere_container();
      if ( sc->needs_rebalance() )
        sc->rebalance();
      
      // debug
      if (0)
      {
        std::cout << "search " << i << " found ";
        _searches[i+1]->sphere_container()->print();
      }
      
      assert( old_dist = _distances[i] ); // yes, assignment inside assert, just to keep it only running when asserts are active
    }
  }

  void add_sphere( size_t si )
  {
    for ( size_t i = 0; i < _searches.size() - 1; ++i )
    {
      _searches[i]->add_sphere(si);
    }
  }
  
  void report_stats(std::ostream &out = std::cout)
  {
    out << "Nested Searches performance: " << std::endl;
    assert( _searches.size() == _distances.size() + 1);
    for ( size_t i = 0; i < _searches.size(); ++i )
    {
      out << i << ": ";
      if (i>0) out << "subset at distance = " << _distances[i-1] << ", ";
      _searches[i]->_stats.report_stats(out);
      Rebalance_Stats *bs = _searches[i]->rebalance_stats();
      if (bs->cpu_time > 0.)
      {
        out << "  Rebalanced " << bs->num_rebalances << " times, in " << bs->cpu_time << " seconds." << std::endl;
      }
    }
  }

};

#endif