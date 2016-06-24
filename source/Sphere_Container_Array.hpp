// Sphere_Container_Array.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// indices into the "Spheres", stored in an array

#ifndef SPHERE_CONTAINER_ARRAY_HPP
#define SPHERE_CONTAINER_ARRAY_HPP

#include <algorithm>
#include "Sphere_Container.hpp"

class Sphere_Container_Array : public Sphere_Container
{
public:
  // big enough to hold small_max_size of the spheres
  Sphere_Container_Array(Spheres &spheres, size_t small_max_size) : Sphere_Container(spheres), _array(0),
    _max_size(small_max_size), _size(0) { new_memory(); }
  // big enough to hold all of the spheres
  Sphere_Container_Array(Spheres &spheres) : Sphere_Container(spheres), _array(0),
    _max_size(spheres.max_num_spheres()), _size(0) { new_memory(); }
  virtual ~Sphere_Container_Array() { delete_memory(); }

  // The first makes the container just big enough for the current items
  Sphere_Container_Array( const Sphere_Container_Array & copy_sca )
  : Sphere_Container( copy_sca._spheres ), _array(0), _max_size( copy_sca.size() ), _size(0)
  {
    new_memory();
    *this = copy_sca; // see copy operator below
  }
  Sphere_Container_Array( const Sphere_Container_Array & copy, size_t max_size );

  // rvalue access only through operator []
  // size_t node = sphere_container_array[9];
  const size_t& operator[](size_t i) const
  {
    assert(i < size());
    return _array[i];
  };
  // lvalue 
  size_t & operator[](size_t i)
  {
    return _array[i];      /* actual access, e.g. return mVector[idx]; */
  };

  // Sphere_Container interface functions
  
  // clear
  // no need to change _array contents, just change size to zero
  virtual void clear()
  { _size = 0;}

  // set the size to be the max size: the entries are not necessarily any good values, however
  // this is useful for using the array as a map
  void enable_all_indices()
  { _size = _max_size; }

  virtual void add_sphere(size_t isphere) 
  { 
    assert(_size < _max_size);
    assert( isphere < _spheres.max_num_spheres() );
    _array[_size++] = isphere; 
  }

  virtual size_t size() const {return _size;}
  size_t max_size() const {return _max_size;}

  Sphere_Container_Array & operator= (const Sphere_Container_Array & other)
  {
    if (this != &other) // protect against invalid self-assignment
    {
      // Modified from comments below: does not allocate new memory

      // 1: allocate new memory and copy the elements
      // int * new_array = new int[other.count];
      assert( _max_size >= other.size() );
      std::copy(other._array, other._array + other.size(), _array);
 
      // 2: deallocate old memory
      // delete [] array;
 
      // 3: assign the new memory to the object
      //array = new_array;
      _size = other.size();
    }
    // by convention, always return *this
    return *this;
  }

  // sort array by ascending index, useful for verification
  void sort()
  {
    std::sort( _array, _array + size() );
  }

  friend class KD_Tree;
  friend class Grid_Container;
  friend class Range_Tree;
protected:
  
  void new_memory()
  {
    delete_memory();
    _array = new size_t[_max_size];
  }
  void delete_memory()
  {
    delete [] _array;
    _array=0;
    _size=0;
  }

  // length of array memory
  size_t _max_size;
  size_t _size;

  size_t *_array;

  virtual size_t next_it( size_t &it ) const
  {
    return ( (it < size()) ? _array[it++] : _spheres.bad_sphere_index() );
  }

};

#endif