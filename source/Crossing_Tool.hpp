// Crossing_Tool.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef CROSSING_TOOL
#define CROSSING_TOOL

#include <algorithm>

#include "Crossing.hpp"
#include "Interval_Tool.hpp"
// see Crossing.hpp for the definition of a Crossing

#include "Spheres.hpp"
#include "Geometry_Tool.hpp"

class Sphere_Container;

// geometry is abstracted, only knows about theta, doesn't know about c, u, v, 
class Crossing_Tool : public Interval_Tool
{
public:
  // convert the piercings into uncovered segments
  // segments shorter than min_length are discarded
  // return true if there are none
  bool uncovered_arcs(Crossings &crossings, const size_t &depth_at_u, std::vector<double> &arcs,
    double &sum_lengths, double min_length = 0.);
protected:
  void sort_crossings( Crossings &crossings )
  { std::sort( crossings.begin(), crossings.end(), theta_less ); }

};

// knows about the domain, etc
class Crossings_Tool : public Crossing_Tool
{
public:
  // Warning: for ghosts, this does not work unless ghosts have been created locally, say by "Search_Ghost::all_near_spheres"
  // or by using permanent ghosts
  Crossings_Tool(Spheres &spheres) :
  _spheres(spheres), _gt( spheres.num_dim() ) {}

  // gather the crossings of the circle of sphere[ci] in the plane of u,v, with the sphere[si]
  // organize by quandrant, in terms of the basis vectors u and v
  void gather_crossings(const double *c, const double rc, double *u, double *v,
                        Sphere_Container *trim_spheres, const double radius_factor,
      Crossings &crossings, size_t &depth);

  // cross_domain not implemented

private:
  const Spheres &_spheres;
  Geometry_Tool _gt;

private:
};

#endif
