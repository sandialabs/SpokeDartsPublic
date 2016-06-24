// Geometry_Tool.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Geometry_Tool.hpp"

#include <cmath>


bool Geometry_Tool::points_of_arc_piercing_sphere(  const double *c, const double rc, 
                                     const double *u, const double *v,
                                     const double *s, double rs,
                                     Crossing &cross_1, Crossing &cross_2, bool &contains_u )
{

  // points
  //
  // vector from c to s
  // double* c_s = new double[_num_dim];
  //
  // point projection of s onto the plane defined by u and v
  // double* sp = new double[_num_dim];
  //
  // vector from c to sp
  // double *c_sp = new double[_num_dim];
  //
  // vector from sp to s
  // double *sp_s = new double[_num_dim];

  // vector from c to sp, normalized to unit length
  // double *c_sp_unit = new double[_num_dim];
  // vector perpendicular to c_sp, in the plane of u and v
  // double *c_sp_perp = new double[_num_dim];
  // point circle-weighted midpoint of edge from c to sp
  //  double *m = new double[_num_dim];
  // vector from c to circle-weighted midpoint of edge from c to sp
  // double *c_m = new double[_num_dim];
  // m_c / r
  // double *c_m_A = new double[_num_dim];
  // vector from c, to intersection of the two circles, the two spheres in the plane of u v
  // normalized to unit distance, rather than distance A, by dividing by A
  // double *i1 = new double[_num_dim];
  // double *i2 = new double[_num_dim];
  
  Geometry_Tool_Workspace *ws = &geometry_tool_workspace;
  double *c_s = ws->c_s;
  double *c_sp = ws->c_sp;
  double *sp_s = ws->sp_s;
  
  // project s to plane of uv, point sp
  // cs = s - c
  subtract( c_s, s, c );
  // vector from c to projection of s
  // c_sp = (c_s * u) u + (c_s * v) v
  const double a_c_sp = dot_product(c_s,u);
  const double b_c_sp = dot_product(c_s,v);
  axpby( c_sp, a_c_sp, u, b_c_sp, v);
  // sp += c_sp + c
  // point_add( sp, c_sp, c);
  
  // d = distance between sp and c
  // const double d_squared = norm_squared( c_sp );
  const double d_squared = a_c_sp * a_c_sp + b_c_sp * b_c_sp; // == norm_squared( c_sp );
  const double d = sqrt(d_squared);
  
  // find radius of the projected sphere, sp_r
  // height of s over sp
  subtract( sp_s, c_sp, c_s);
  const double r_sp_squared = rs * rs - norm_squared( sp_s );
  const double r_sp = sqrt( r_sp_squared );
  
  // if the circles intersect
  //   first, the projected disk close enough to the center
  //   second, the disk is not so close it is inside the rc-radius circle.
  if ( d < r_sp + rc && r_sp + d > rc )
  {
    
    // hi = distance from sp_c line to each intersection point of the circles
    // see http://mathworld.wolfram.com/Circle-CircleIntersection.html formula for "a"
    assert( (-d + r_sp - rc) * (-d - r_sp + rc) * (-d + r_sp + rc) * ( d + r_sp + rc) > 0 );
    const double hi = sqrt( (-d + r_sp - rc) * (-d - r_sp + rc) * (-d + r_sp + rc) * ( d + r_sp + rc) ) / (2. * d);
    assert( hi < 1.0001 * r_sp );
    // fraction of r the intersection point is above m
    const double hi_rc = hi / rc;
    assert( hi_rc < 1.0001 );
    
    // m: circle-weighted midpoint of edge from c to sp
    // wolfram's formula for x
    const double c_m_distance = (d_squared - r_sp_squared + rc*rc) / (2. * d);
    // m = (c_m_distance / d) * c_sp + c, where d = norm(c_sp)
    // point_axpy( m, c_m_distance / d, c_sp, c);
    // m_c = (c_m_distance / d) * c_sp, where d = norm(c_sp)
    // point_multiply( m_c, c_m_distance / d, c_sp );
    // point_multiply( c_sp_unit, c_sp, 1. / d );
    // assert( fabs(norm_squared(c_sp_unit) - 1.) < 0.001 );
    // m_frac = distance from c to m divided by rc
    const double m_frac =  c_m_distance / rc;
    assert( m_frac <= 1.0001 );
    // point_multiply( c_m_A, c_sp_unit, m_frac );
    
    // (a_c_sp/d, b_c_sp/d) is normalized to unit length
    const double a_c_m_unit = a_c_sp / d;
    const double a_c_m = m_frac * a_c_m_unit;
    const double b_c_m_unit = b_c_sp / d;
    const double b_c_m = m_frac * b_c_m_unit;
    // a^2 + b^2 = m_frac^2
    assert( fabs(a_c_m_unit * a_c_m_unit + b_c_m_unit * b_c_m_unit - 1.) < 0.001);
    assert( fabs(a_c_m * a_c_m + b_c_m * b_c_m - m_frac * m_frac) < 0.001);
    
    //  c_sp_perp
    // perp of (a u + b v ) is ( -b u + a v )
    // point_axpby( c_sp_perp, -b_c_sp, u, a_c_sp, v);
    
    // i1 and i2 have been normalized to unit radius
    // i1 = m_c + hi * c_sp_perp = a_c_m u + b_c_m v + (hi * -b_c_sp) u  + (hi *  a_c_sp) v
    // i2 = m_c - hi * c_sp_perp = a_c_m u + b_c_m v + (hi *  b_c_sp) u  + (hi * -a_c_sp) v
    //
    // point_axpy( i1,  hi / r, c_sp_perp, m_c_r);
    // point_axpy( i2, -hi / r, c_sp_perp, m_c_r);
    
    // find the (a,b) pairs of i1 and i2, in the basis u and v, of unit length, not of rc length
    // use the formula for perpendicular vectors for the second column below
    const double a_i1 = a_c_m - hi_rc * b_c_m_unit;
    const double b_i1 = b_c_m + hi_rc * a_c_m_unit;
    const double a_i2 = a_c_m + hi_rc * b_c_m_unit;
    const double b_i2 = b_c_m - hi_rc * a_c_m_unit;
    assert( a_i1 > -1.01 );
    assert( a_i1 <  1.01 );
    assert( b_i1 > -1.01 );
    assert( b_i1 <  1.01 );
    assert( a_i2 > -1.01 );
    assert( a_i2 <  1.01 );
    assert( b_i2 > -1.01 );
    assert( b_i2 <  1.01 );
    assert( fabs( a_i1 * a_i1 + b_i1 * b_i1 - 1. ) < 0.001 );
    assert( fabs( a_i2 * a_i2 + b_i2 * b_i2 - 1. ) < 0.001 );
    
    // angle from u, in [0,2pi]
    double theta_1 = atan2(b_i1,a_i1); // slow? fast proxy?
    double theta_2 = atan2(b_i2,a_i2); // slow? fast proxy?
    if (theta_1 < 0.)
      theta_1 += 2. * PI;
    if (theta_2 < 0.)
      theta_2 += 2. * PI;
    assert( theta_1 >= 0. );
    assert( theta_1 < 2. * PI );
    assert( theta_2 >= 0. );
    assert( theta_2 < 2. * PI );
    
    // order so theta_1 < theta_2p
    bool forward1 = false;
    contains_u = false;
    const double theta_2p = (theta_1 < theta_2) ? theta_2 : theta_2 + 2. * PI;
    assert(theta_2p > theta_1);
    if (theta_2p - theta_1 < PI)
    {
      assert( (theta_1 < theta_2) || (theta_1 - theta_2 >= PI) );
      forward1 = true;
      if ( !(theta_1 < theta_2) )
        contains_u = true;
    }
    else
    {
      forward1 = false;
      if (theta_1 < theta_2)
        contains_u = true;
    }
    if (forward1)
    {
      cross_1 = Crossing(a_i1, b_i1, theta_1,  forward1);
      cross_2 = Crossing(a_i2, b_i2, theta_2, !forward1);
    }
    else
    {
      cross_2 = Crossing(a_i1, b_i1, theta_1,  forward1);
      cross_1 = Crossing(a_i2, b_i2, theta_2, !forward1);      
    }
    return true;
  }
  // else the circles didn't intersect and there is nothing to do
  return false;
}


bool Geometry_Tool::points_of_line_piercing_sphere( const double *c, const double *u, const double *s,
                                     double &x, double &y, const double radius_factor )
{
  Geometry_Tool_Workspace *ws = &geometry_tool_workspace;
  double *c_s = ws->c_s;

  const double r = s[num_dim()] * radius_factor;

  // c_s, vector from c to s
  subtract( c_s, s, c );

  // distance from c to (projecction of s onto the line). could be negative
  const double a_c_sp = dot_product(c_s,u);

  // distance from (projection of s onto line) to crossing points,
  const double pm_squared = a_c_sp * a_c_sp + r*r - norm_squared( c_s );

  // if sphere radius is too small to intersect the ray, quit
  if (pm_squared <= 0.)
    return false;
  
  const double pm = sqrt( pm_squared );
  
  x = a_c_sp - pm;
  y = a_c_sp + pm;
  return true;
}


void Geometry_Tool::separating_hyperplane( const double *c, const double *s, double *p, double *n, double &pdist )
{
  // p = c + s / 2
  assign(p,c);
  add(p,s);
  multiply(p, 0.5);
  
  // n = s - c, normalized
  assign(n, s);
  subtract(n, c);
  pdist = normalize(n) * 0.5;
}

