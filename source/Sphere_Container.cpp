// Sphere_Container.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Sphere_Container.hpp"
#include "Const_Sphere_Container_Iterator.hpp"

#include <stdio.h>

void Sphere_Container::print( std::string name, std::ostream &out ) const
{
  print_distances(0, name, out);
}


void Sphere_Container::print_distances(const double *c, std::string name, std::ostream &out, size_t max_lines ) const
{
  Const_Sphere_Container_Iterator iiter( this );
  out << name << " sphere container " << this << " contains " << size() << " spheres:" << std::endl;
  if (c)
  {
    out << " distances from point ";
    _spheres._pt.print_sphere(c,out);
    out << std::endl;
  }
  size_t num_lines(0);
  for (size_t i = iiter.first(); i != iiter.bad_sphere_index(); i = iiter.next())
  {
    out << "sphere " << i << ": ";
    _spheres._pt.print_sphere( _spheres[i], out );
    if (c)
    {
      const double dist = sqrt( _spheres._pt.distance_squared(c, _spheres[i]) );
      out << " dist:" << dist;
    }
    out << std::endl;
    // avoid flooding the screen with too much data
    if (max_lines && (++num_lines > max_lines))
        break;
  }
}


bool Sphere_Container::verify_min_distances( std::ostream &out ) const
{
  if (!size())
    return true;
  bool all_distances_big = true;
  double r, r_tol, r_squared, r_squared_tol;
  get_min_distances( r, r_tol, r_squared, r_squared_tol );

  Const_Sphere_Container_Iterator iiter( this );
  for (size_t i = iiter.first(); i != iiter.bad_sphere_index(); i = iiter.next())
  { 
    const double *p = _spheres[i];
    // start counting j after i, to avoid 1/2 the computation
    Const_Sphere_Container_Iterator jiter( iiter );
    for (size_t j = jiter.next(); j != jiter.bad_sphere_index(); j = jiter.next())
    { 
      if ( !verify_min_distances_loc( p, _spheres[j], i, j,
            r, r_tol, r_squared, r_squared_tol, out) )
        all_distances_big = false;
    }
  }
  return all_distances_big;
}

bool Sphere_Container::verify_min_distances( const double *p, const size_t pi, std::ostream &out ) const
{
  if (!size())
    return true;
  bool all_distances_big = true;
  double r, r_tol, r_squared, r_squared_tol;
  get_min_distances( r, r_tol, r_squared, r_squared_tol );

  Const_Sphere_Container_Iterator jiter( this );
  for (size_t j = jiter.first(); j != jiter.bad_sphere_index(); j = jiter.next())
  { 
    if ( !verify_min_distances_loc( p, _spheres[j],  pi, j,
            r, r_tol, r_squared, r_squared_tol, out) )
      all_distances_big = false;
  }
  return all_distances_big;
}

bool Sphere_Container::verify_min_distances( const double *p, const double *q,
                                            const size_t pi,const size_t qi,
                                            std::ostream &out ) const
{
  double r, r_tol, r_squared, r_squared_tol;
  get_min_distances( r, r_tol, r_squared, r_squared_tol );
  return verify_min_distances_loc( p, q,
                                  pi, qi,
                                  r, r_tol, r_squared, r_squared_tol,
                                  out);
}

bool Sphere_Container::verify_min_distances_loc( const double *p, const double *q, 
                                                size_t pi, size_t qi,
                                                double &r, double&r_tol, double &r_squared, double &r_squared_tol,
                                                std::ostream &out ) const
{
  // ? check if p or q is null? Is that valid or not?
  if ( p != q )
  {
    const double distance_squared = _spheres._pt.distance_squared( p, q );
    if ( distance_squared < r_squared_tol )
    {
      const double distance = sqrt( distance_squared );
      out << "Conflict distance error,";
      out << " disks at distance:" << distance << " < radius:" << r << " (threshold " << r_tol << ")" << std::endl;
      if ( pi != _spheres.bad_sphere_index() )
        out << " i:" << pi;
      _spheres._pt.print_point( p, out); out << std::endl;
      if (qi != _spheres.bad_sphere_index() )
        out << " j:" << qi;
      _spheres._pt.print_point( q, out); out << std::endl;
      return false;
    }
  }
  return true;
}

void Sphere_Container::get_min_distances( double &r, double&r_tol, double &r_squared, double &r_squared_tol) const
{
  // for now, assume all spheres have the same radii
  r = _spheres.radius(0);
  r_tol = r * (1. - 1.e-4);
  r_squared = r*r;
  r_squared_tol = r_tol * r_tol;
}
