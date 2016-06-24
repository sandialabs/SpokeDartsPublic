// Crossing.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// a crossing is the point where a cirular arc crosses a sphere, along with whether the sphere is entered or left in
// the forward direction of the circle.
// This defines intervals of arcs inside the sphere.
// These piercings can be used to efficiently determine the non-covered or covered arcs of the circle

// compare Crossings for arcs of circles to Piercings for segments of lines

#ifndef CROSSING
#define CROSSING

#include <vector>

class Crossing
{
public:
  // the crossing point is r(a*u + b*v) + c
  // it makes angle theta with the vector u
  double a, b, theta;
  bool forward; // true if this is the beginning of a covered arc, false if it is the end
  
  Crossing(): a(0.), b(0.), theta(0.), forward(true) {;}
  Crossing(const double a_in, const double b_in, const double theta_in, const bool forward_in)
  : a(a_in), b(b_in), theta(theta_in), forward(forward_in) {}
  
  // operators for sorting by angle
  bool operator() (Crossing  i, Crossing  j) { return ( i.theta <  j.theta);}
  bool operator() (Crossing *i, Crossing *j) { return (i->theta < j->theta);}
};

// pointer to Crossing might be better
typedef std::vector< Crossing > Crossings;
// e.g. 
// Crossings crossings;
// crossings.push_back( Crossing(a_i1, b_i1, theta_1,  forward1) );

// for sorting by angle
// the xcode compiler gives a spurious warning that this function is unused, but it actually is and won't compile if it is commented out.
static bool theta_less(Crossing  i, Crossing  j) { return ( i.theta <  j.theta); }
// e.g.  std::sort( crossings.begin(), crossings.end(), theta_less );


#endif