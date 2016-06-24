// Sphere_Container_Iterator.hpp 

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// Abstract class
//
// Interface to an simple iterator over a generic container
// Should not change the container in any way, const iterator conceptually,

#ifndef SPHERE_CONTAINER_ITERATOR_HPP
#define SPHERE_CONTAINER_ITERATOR_HPP

#include "Spheres.hpp"
class Sphere_Container;

class Sphere_Container_Iterator
{
public:  
  Sphere_Container_Iterator(Spheres &spheres) : _spheres(spheres), _it(0) {}
  virtual ~Sphere_Container_Iterator() {}

  // interface:
  //
  // reset_iterator()
  // first()
  // next()
  // bad_sphere_index()
  //
  // idiom:
  // Sphere_Container_Iterator *sci = Sphere_Container, or Const_Sphere_Container_Interator
  //  for (size_t i = sci->first(); i != sci->bad_sphere_index(); i = sci->next())
  //  { _spheres[i]; }
  // or:  double * s = _spheres[ next() ]; 
  //
  // first and next return the current value of it, and it is advanced to the next valid value or to bad_sphere_index

  // non-virtual
  // these interfaces just call the const version with its own index
  void reset_iterator() { reset_iterator_it(_it); }

  // after the last, next returns _spheres.bad_sphere_index()
  size_t next() { return next_it(_it); }
  
  size_t first()
  {
    reset_iterator_it(_it);
    return next_it(_it);
  }

  size_t it() const {return _it;}
  
  size_t bad_sphere_index() const { return _spheres.bad_sphere_index(); }

  size_t num_dim() const { return _spheres.num_dim(); }

  Spheres &_spheres; // for _spheres[ next() ];
  Spheres &spheres() {return _spheres;}
  const Spheres &spheres() const {return _spheres;}
  
  friend class Const_Sphere_Container_Iterator;

protected:
  // used for defining an implementation-specific version of the above that is const, 
  // for a const-like integer-based iterator

  // next must be defined by derived classes, and may be non-trivial
  // reset_iterator might need to be redefined if the first container index is not 0
  virtual void reset_iterator_it( size_t &it ) const { it = 0; }
  virtual size_t next_it( size_t &it ) const = 0;
  // first() is not virtual, it doesn't need to be defined again
  
  size_t _it;
};

#endif