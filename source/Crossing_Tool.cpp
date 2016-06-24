// Crossing_Tool.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Crossing_Tool.hpp"
#include "Sphere_Container.hpp"

void Crossings_Tool::gather_crossings(const double *c, const double rc, double *u, double *v,
                                      Sphere_Container *trim_spheres, const double radius_factor,
      Crossings &crossings, size_t &depth)
{
  bool contains_u(false);
  depth = 0;
  Crossing cross_1, cross_2;
  for (size_t i = trim_spheres->first(); i != trim_spheres->bad_sphere_index(); i = trim_spheres->next())
  {
    const double *s = _spheres[i];
    double rs = _gt.radius(s) * radius_factor;
    _gt.points_of_arc_piercing_sphere( c, rc, u, v, s, rs, cross_1, cross_2, contains_u ); 
    crossings.push_back( cross_1 );
    crossings.push_back( cross_2 );
    if (contains_u)
      depth++;
  }
}

bool Crossing_Tool::uncovered_arcs(Crossings &crossings, const size_t &depth_at_u, std::vector<double> &arcs,
    double &sum_lengths, double min_length)
{

  sort_crossings( crossings );
  
  // traverse to find uncovered segments
  size_t current_depth = depth_at_u;
  // first interval containing u will be covered, as otherwise we would have used the spoke along it
  assert(depth_at_u != 0);
  sum_lengths = 0;
  for ( size_t i = 0; i < crossings.size(); ++i )
  {
    if (crossings[i].forward)
    {
      ++current_depth;
      // did I just finish an interval?
      if ( current_depth == 1)
      {
        const double interval_length = crossings[i].theta - crossings[i-1].theta;
        if (interval_length > min_length)
        {
          sum_lengths += interval_length;
          arcs.push_back(crossings[i-1].theta);
          arcs.push_back(crossings[i].theta);
        }
      }
    }
    else
    {
      assert( current_depth > 0 );
      --current_depth;
    }
  }
  // final depth should be the same as the starting depth, since its a circle
  assert(current_depth == depth_at_u);
  // add the segment crossing u if it is uncovered
  if (current_depth == 0)
  {
    arcs.push_back(crossings.back().theta);
    arcs.push_back(crossings.front().theta + 2 * PI);
  }
  return !arcs.empty();
}
