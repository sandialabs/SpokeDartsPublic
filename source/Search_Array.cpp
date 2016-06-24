// Search_Array.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include <limits>

#include "Search_Array.hpp"
#include "Sphere_Container.hpp"
#include "Sphere_Container_Array.hpp"
#include "Global_Container.hpp"
#include "Ghost_Global_Container.hpp"

Search_Array::Search_Array(Spheres &spheres, bool set_global ) : _sphere_container(0), Search_Structure(spheres.num_dim()), _this_owns_container(true)
{
  _is_global = set_global;
  if (is_global())
    _sphere_container = Ghost_Global_Container::new_global_container(spheres);
  else
  {
    _sphere_container = new Sphere_Container_Array(spheres);
    _this_owns_container = true;
  }
}

void Search_Array::change_container( Sphere_Container *new_container, bool set_global )
{
  Spheres *spheres = & sphere_container()->_spheres;
  if (_this_owns_container)
    delete _sphere_container;
  if (new_container)
  {
    _sphere_container = new_container;
    _this_owns_container = false;
  }
  else
  {
    _is_global = set_global;
    if (is_global())
      _sphere_container = Ghost_Global_Container::new_global_container(*spheres);
    else
      _sphere_container = new Sphere_Container_Array(*spheres);
    
    _this_owns_container = true;
  }
}

double Search_Array::nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold )
{
  _stats.start_clock();
  double thresh_squared = ( distance_threshold < sqrt( std::numeric_limits<double>::max() ) ) ? 
    distance_threshold * distance_threshold : std::numeric_limits<double>::max();
  near_sphere = Spheres::bad_sphere_index();

  for (size_t i = sphere_container()->first(); i != sphere_container()->bad_sphere_index(); i = sphere_container()->next())
  {
    const double *q = spheres()[i];
    if (q != p)
    {
      const double pair_distance = _tt.distance_squared( p, q);
      if (thresh_squared >= pair_distance )
      {
        thresh_squared = pair_distance;
        near_sphere = i;
      }
    }
  }

  _stats.collect_stats(sphere_container()->size(), sphere_container()->size(), 1);
  return ( sqrt( thresh_squared ) );
}


bool Search_Array::no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor, const double dist)
{
  _stats.start_clock();
   // tolerance will skip some almost-close spheres
  const double thresh = dist + spheres().radius(0) * radius_factor -1E-10;
  const double thresh_squared = thresh * thresh;
  size_t v = 0;
  for (size_t i = sphere_container()->first(); i != sphere_container()->bad_sphere_index(); i = sphere_container()->next())
  {
    ++v;
    const double pair_distance = _tt.distance_squared( p, spheres()[i]);
    if (pair_distance < thresh_squared)
    {
      near_sphere = i;
      _stats.collect_stats(sphere_container()->size(), v, 1);
      return false; 
    }
  }
  _stats.collect_stats(sphere_container()->size(), v, 0);
  return true;
}

void Search_Array::all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius )
{
  _stats.start_clock();
  sub_container->clear();
  double thresh = dist;
  if ( subtract_radius )
    thresh += spheres().radius(0);
  const double thresh_squared = thresh * thresh;
  for (size_t i = sphere_container()->first(); i != sphere_container()->bad_sphere_index(); i = sphere_container()->next())
  {
    if ( spheres()[i] != p )
    {
      const double pair_distance = _tt.distance_squared( p, spheres()[i]);
      if (pair_distance < thresh_squared)
      {
        sub_container->add_sphere(i);
      }  
    }
  }
  _stats.collect_stats(sphere_container()->size(), sphere_container()->size(), sub_container->size());
}

void Search_Array::trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2, const double radius_factor  )
{
  _stats.start_clock();
  size_t true_trims(0);
  for (size_t i = sphere_container()->first(); i != sphere_container()->bad_sphere_index(); i = sphere_container()->next())
  {
    if (_tt.anchored_trim( c, u, spheres()[i], A, A_1, A_2, radius_factor ) )
      ++true_trims;
    assert( A_1 <= A_2 ); // but either could be negative
    assert( fabs(A_1) >= r * (0.98) - 1e-6); // the segment must be outside the initial sphere
    assert( fabs(A_2) >= r * (0.98) - 1e-6); 
  }
  _stats.collect_stats(sphere_container()->size(), sphere_container()->size(), true_trims);
}
