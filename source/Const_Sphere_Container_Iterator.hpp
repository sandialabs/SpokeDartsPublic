// Const_Sphere_Container_Iterator.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// Basic const iterator, over an existing iterator of a container

#ifndef CONST_SPHERE_CONTAINER_ITERATOR_HPP
#define CONST_SPHERE_CONTAINER_ITERATOR_HPP

#include "Sphere_Container_Iterator.hpp"

class Const_Sphere_Container_Iterator : public Sphere_Container_Iterator
{
public:
  // a kind of copy constructor
  // note the current position of the iterator is also copied.
  // use first() or reset_iterator() if you want to start from the beginning
  Const_Sphere_Container_Iterator(const Sphere_Container_Iterator *sci) : Sphere_Container_Iterator(sci->_spheres), _sci{sci} { _it = sci->it(); }
  virtual ~Const_Sphere_Container_Iterator() {}

protected:

  const Sphere_Container_Iterator *_sci;

  virtual void reset_iterator_it( size_t &it ) const
  {
    _sci->reset_iterator_it( it );
  }

  virtual size_t next_it( size_t &it ) const
  {
    return _sci->next_it( it );
  }

};

#endif