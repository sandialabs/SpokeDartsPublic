// Trimming.hpp
// geometric tool for taking a spoke and reducing its extent by other spheres or the domain boundary

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef TRIMMING_TOOL
#define TRIMMING_TOOL

#include "Geometry_Tool.hpp"
class Domain;

// All the geometry operations for trimming by one sphere.
// Not the tree traversal for trimming 

class Trimming_Tool : public Geometry_Tool
{
  public:

  Trimming_Tool(size_t num_dim) : Geometry_Tool(num_dim) {}

  // for line from c in direction u for length t,
  // reduce t so that the line does not leave the domain boundary
  // return true if t was actually reduced
  bool trim_by_domain( const double *c, const double *u, double &t, const double *xmax, const double *xmin );
  bool trim_by_domain( const double *c, const double *u, double &t, Domain *domain );

//zzyk todo: make axis-aligned line versions of the above, for kd-darts

  // spoke line goes from c in direction u between A_1 and A_2, with anchor A between them
  // s is the trimming sphere, dialated by the radius_factor
  // return true if the spoke was actually reduced
  bool anchored_trim(const double *c, const double *u, const double *s, double A, double &A_1, double &A_2, const double radius_factor );

  // Trim by the Voronoi hyperplane separating the two spheres.
  // Retain only the portion of the spoke closer to c than s.
  // Assume equal weights for now.
  // Caller can tell if the anchor was trimmed:
  //   bool anchor_was_trimmed = (A_2 < A) || (A_1 > A);
  // This runs faster if workspace p and n are passed in
  bool hyperplane_trim(const double *c, const double *u, const double *s, double A, double &A_1, double &A_2,
                       double *p=0, double *n=0);


private:

};

#endif