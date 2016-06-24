// Trimming_Tool.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Trimming_Tool.hpp"
#include "Domain.hpp"

bool Trimming_Tool::trim_by_domain( const double *c, const double *u, double &t, Domain *domain)
{
  if (domain->is_periodic())
    return false;
  return trim_by_domain( c, u, t, domain->xmax(), domain->xmin() );
}


bool Trimming_Tool::trim_by_domain( const double *c, const double *u, double &t, const double *xmax, const double *xmin )
{
  const double t_init = t;
  for (size_t idim = 0; idim < num_dim(); idim++)
  {
    const double x = c[idim] + t * u[idim];
    
    // t will be the minimum of all these quantities over all idim, successive reduction
    if      (x > xmax[idim]) t = (xmax[idim] - c[idim]) / (u[idim] * (1. + 1e-10));
    else if (x < xmin[idim]) t = (xmin[idim] - c[idim]) / (u[idim] * (1. + 1e-10));
  }
  // fix roundoff
  if ( fabs(t) > fabs(t_init) ) 
    t = t_init;

  // debug
  // check that domain trimming worked
  /*
  for (size_t idim = 0; idim < _num_dim; idim++)
  {
    double x = c[idim] + t_abs * u[idim];
    assert (x <= xmax[idim]);
    assert (x >= xmin[idim]);
  }
   */

  return (t < t_init);
} 

bool Trimming_Tool::anchored_trim(const double *c, const double *u, const double *s, double A, double &A_1, double &A_2, const double radius_factor )
// points or vectors
// p -> c  // center of sphere
// q ->    // dart point
//      u  // unit vector from c towards a, the spoke direcion
//      s  // sphere we use to trim this line, dilated by radius_factor
{
  double x, y;
  bool ret_val(false);
  if ( points_of_line_piercing_sphere( c, u, s, x, y, radius_factor) )
  {
    assert( y >= x );

    // disc is left
    if ( y < A )
    {
      if (y > A_1)
      {
        A_1 = y;
      }
    }
    // disc is right
    else if (x > A)
    {
      if (x < A_2)
        A_2 = x;
    }
    // disk straddles A, by roundoff
    // then get rid of one side or the other of the interval
    else
    {
      if (A - x < y - A)
        A_2 = A;
      else
        A_1 = A;
    }
    assert(A_1 <= A);
    assert(A_2 >= A);
    ret_val = true;
  }
  
  // debug
  // check that the points of the trimmed segment are actually outside the sphere s
  if (_verification_level >= 1 )
  {
    double *a = new_point();
    
    // point at A
    axpy( a, A, u, c );
    assert( distance_squared(s, a) > radius(s) * radius(s) * 0.9999 );
    
    // point at A_1
    axpy( a, A_1, u, c );
    assert( distance_squared(s, a) > radius(s) * radius(s) * 0.9999 );
    
    // point at A_2
    axpy( a, A_2, u, c );
    assert( distance_squared(s, a) > radius(s) * radius(s) * 0.9999 );
    
    delete_point(a);
  }
  
  return ret_val;
}


bool Trimming_Tool::hyperplane_trim(const double *c, const double *u, const double *s, double A, double &A_1, double &A_2, double *p, double *n)
{
  // point and normal defining hyperplane
  const bool own_points = (p==0 || n==0);
  if (own_points)
  {
    p = new_point();
    n = new_point();
  }
  double mid_dist(0.);
  separating_hyperplane( c, s, p, n, mid_dist );
  mid_dist /= radius(c);
    
  const double dot = dot_product(u,n);
  
  // does intersection take place far away?
  if ( mid_dist > fabs(A_2 * dot) )
    return false;
  
  // the above check ensures we don't divide by zero here
  const double Anew = mid_dist / dot;
  
  bool was_trimmed = false;
  if (Anew > 0)
  {
    if ( Anew < A_2 )
    {
      was_trimmed = true;
      A_2 = Anew;
    }
  }
  else /*Anew < 0*/
  {
    if ( A_1 < Anew )
    {
      was_trimmed = true;
      A_1 = Anew;
    }
  }

  // caller can tell if the anchor was trimmed:
  // bool anchor_was_trimmed =  ( A_2 < A ) || (A_1 > A);
  
  if (own_points)
  {
    delete_point(p);
    delete_point(n);
  }
  
  return was_trimmed;
}
