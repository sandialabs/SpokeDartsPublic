// Search_Array.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// The search routines just iterate over all the points
// it has a container, in this case a simple array or list-like structure

#ifndef SEARCH_ARRAY_HPP
#define SEARCH_ARRAY_HPP

#include "Trimming_Tool.hpp"
#include "Search_Structure.hpp"
#include "Spheres.hpp"
#include "Sphere_Container.hpp"

// This is the geometric search methods performed over a simple linear array (vector)

class Search_Array : public Search_Structure
{
public:

  // creates a global container for the spheres, and searches over these
  Search_Array(Spheres &spheres, bool set_global = true);
  // uses the passed in container, and searches the spheres in that container
  Search_Array(Spheres &spheres, Sphere_Container &sphere_container) 
  : _sphere_container(&sphere_container), _this_owns_container(false), Search_Structure(spheres.num_dim()) {}
  // if new_container==NULL, then create and use a global one instead
  void change_container( Sphere_Container *new_container, bool set_global = true );

  // See Search_Structure interface class for descriptions of these functions
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() );

  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.);

  virtual void all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false );

  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
    const double radius_factor = 1. );

//  virtual void trim_circle_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
//    const double radius_factor = 1. );

  virtual Sphere_Container *sphere_container() {return _sphere_container;}
  
  virtual ~Search_Array()
  {
    if (_this_owns_container)
      delete _sphere_container;
    _sphere_container = 0;
  }
  
protected:

private:
  Search_Array();

  Sphere_Container *_sphere_container;
  // if true, then this new'ed the container, and should delete it on destruction
  bool _this_owns_container;

  Spheres &spheres() {return sphere_container()->_spheres;}


};

#endif