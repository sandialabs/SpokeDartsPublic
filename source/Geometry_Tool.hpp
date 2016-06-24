// Geometry_Tool.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef GEOMETRY_TOOL
#define GEOMETRY_TOOL

#include "Point_Tool.hpp"
#include "Crossing.hpp"

class Geometry_Tool : public Point_Tool
{
public:

  // constucting a Geometry_Tool is expensive because of dynamic memory allocation
  Geometry_Tool(size_t num_dim)
  : Point_Tool(num_dim), geometry_tool_workspace(this) {}

  // compute the crossings points of a circle with a sphere
  // the circle is centered at c with radius rc, in the plane of orthogonal unit vectors u,v,
  // the sphere is centered at s with radius rs
  // The two points are returned in the Crossings, ordered so that cross_1 is forward
  // crosses_u = if s contains the point c + rc * u  
  // returns whether the sphere actually crossed the circle
  // replaces gather_crossings
  bool points_of_arc_piercing_sphere(  const double *c, const double rc, 
                                       const double *u, const double *v,
                                       const double *s, double rs,
                                       Crossing &cross_1, Crossing &cross_2, bool &contains_u );
  // c and s are spheres, use their natural radii
  bool points_of_arc_piercing_sphere(  const double *c,  
                                       const double *u, const double *v,
                                       const double *s, 
                                       Crossing &cross_1, Crossing &cross_2, bool &crosses_u )
  {
    return points_of_arc_piercing_sphere(c, radius(c), u, v, s, radius(s), cross_1, cross_2, crosses_u);
  }
  // not implemented, find the two points where the circle leaves the domain
  //  void gather_crossings_axis( const size_t ci, const double A, const double *u, const double *v,
  //                              double * xmin, double * xmax,
  //                              Crossings &crossings, size_t &depth );

  // return the two points x and y where a line pierces the sphere s
  // the line goes in direction u from point p, 
  // the sphere s is dialated by the radius_factor
  // return the distances x and y from p to where the line pierces the sphere
  // returns false if the line misses the sphere
  bool points_of_line_piercing_sphere( const double *p, const double *u, const double *s, double &x, double &y,
    const double radius_factor = 1.);
  
  // return midpoint p and normal n of the hyperplane separating points c and s
  void separating_hyperplane( const double *c, const double *s, double *p, double *n, double &pdist );

private:
  Geometry_Tool();

  
  
  // data used within methods, placed here to avoid dynamic memory allocation
  class Geometry_Tool_Workspace
  {
  public:
    double *c_s; // point
    double *c_sp; // point
    double *sp_s; // point
    
    ~Geometry_Tool_Workspace()
    {
      if (_pt)
      {
        _pt->delete_point(c_s);
        _pt->delete_point(c_sp);
        _pt->delete_point(sp_s);
      }
    }
    Geometry_Tool_Workspace(Geometry_Tool *gt)
    : _pt( new Point_Tool(gt->num_dim()))
    {
      c_s  = _pt->new_point();
      c_sp = _pt->new_point();
      sp_s = _pt->new_point();
    }
    
    Point_Tool * _pt;
  };
  Geometry_Tool_Workspace geometry_tool_workspace;
  
};





#endif
