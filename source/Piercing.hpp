// Piercing.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// a piercing is the point where a line crosses a sphere, along with whether the sphere is entered or left in
// the forward direction of the line.
// This defines intervals inside the sphere.
// These piercings can be used to efficiently determine the non-covered or covered segments of the line

// compare Crossings for arcs of circles to Piercings for of segments of a line

#ifndef PIERCING
#define PIERCING

#include <algorithm>

#include "Interval_Tool.hpp"
#include "Spheres.hpp"
#include "Trimming_Tool.hpp"
class Sphere_Container;

class Piercing
{
public:
  double p;
  bool forward; // true if this is the beginning of a covered segment, false if it is the end
  
  Piercing(): p(0.), forward(true) {;}
  Piercing(const double p_in, const bool forward_in)
  : p(p_in), forward(forward_in) {}
  
  // operators for sorting by order along on the ray
  bool operator() (Piercing  i, Piercing j) { return ( i.p <  j.p);}
  bool operator() (Piercing *i, Piercing *j) { return (i->p < j->p);}
};

static bool p_less(Piercing  i, Piercing  j) { return ( i.p <  j.p);}

typedef std::vector< Piercing > Piercings;
// Piercings.push_back( Piercing(a, forward) );
// a vector of *pointers* to Piercing might perform better

class Piercing_Tool : public Interval_Tool
{
public:
  // add a piercing (covered segment of a line) to the lists of piercings
  // the line segment of interest goes from distance A to B,
  // the piercing segment is from x to y
  // depth is the number of piercings covering A
  // return true if the entire line is covered by a single sphere
  bool two_piercings(const double A, const double B, const double x, const double y,
                     Piercings &piercings, size_t &depth_at_A);

  // convert the piercings into uncovered segments
  // segments shorter than min_length are discarded
  // return true if there are none
  bool uncovered_segments(Piercings &piercings, const size_t &depth_at_A, double A, double B,
                          std::vector<double> &segments, double &sum_lengths, double min_length = 0.);

protected:
  void sort_piercings( Piercings &piercings )
  { std::sort( piercings.begin(), piercings.end(), p_less ); }
};


// this is the stuff that deals with the set of spheres and the domain
class Piercings_Tool : public Piercing_Tool
{
public:
  // Warning: for ghosts, this does not work unless ghosts have been created locally, say by "Search_Ghost::all_near_spheres"
  // or by using permanent ghosts
  Piercings_Tool( Spheres &spheres ) : _spheres(spheres), _tt( spheres.num_dim() ) {}

  // find all the piercings made by this array
  // return true if the segment is covered, or if quit_if_segment_touched and it contains any piercings.
  // If Piercings *piercings_m is passed in, then the segment is considered double-ended, as the two sides of a line spoke, from -A to -B
  //   and we only return true if both sides would return true
  bool gather_piercings( Sphere_Container *pierced_spheres,
    const double *p, const double *u, const double A, const double B, 
    Piercings *piercings_p, size_t *depth_at_A_p,
    Piercings *piercings_m = 0, size_t *depth_at_A_m = 0,
    bool quit_if_segment_touched = false);

  // add crossings where the line leaves the domain, for non-periodic domains
  bool pierce_domain( 
    const double *xmax, const double *xmin,
    const double *p, const double *u, const double A, const double B,  
    Piercings *piercings_p, size_t *depth_at_A_p,
    Piercings *piercings_m = 0, size_t *depth_at_A_m = 0,
    bool quit_if_segment_touched = false);

private:
  void pierce_domain_oneside(
    const double *xmax, const double *xmin,
    const double *p, const double *u, const double A, const double B, 
    Piercings *piercings, size_t *depth_at_A,
    bool &touched, bool &covered);

  const Spheres &_spheres;
  Trimming_Tool _tt;

};


#endif
