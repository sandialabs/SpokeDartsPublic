// Global_Container.hpp 

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// contains all the real spheres, always, so no need to store indices

#ifndef GLOBAL_CONTAINER
#define GLOBAL_CONTAINER

#include "Sphere_Container.hpp"
#include "Spheres.hpp"

class Global_Container: public Sphere_Container
{
public:
  Global_Container(Spheres &spheres) : Sphere_Container(spheres) {}

  // Sphere_Container interface functions
  // clear() do nothing, the underlying spheres will have to be emptied instead
  virtual void clear() {}

  // add_sphere() do nothing, it will have to be added to the underlying spheres instead
  virtual void add_sphere(size_t isphere) {}

  virtual size_t size() const 
  { return _spheres.num_real_spheres(); }
  

  // set it up so that the next call to next_it returns i, if i is a valid value, otherwise bad_index
  virtual void set( size_t i )
  {
    _it = i;
  }
  
protected: 

  virtual size_t next_it(size_t &it) const
  { return ( (it < _spheres.num_real_spheres()) ? it++ : _spheres.bad_sphere_index() ); }

};

#endif
