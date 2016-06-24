// KD_Tree.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// This is the data for the tree

#ifndef KD_TREE
#define KD_TREE

#include "Sphere_Container.hpp"
#include "Point_Tool.hpp"
#include "Sphere_Container_Array.hpp"

class KD_Tree : public Sphere_Container
{

public:

  KD_Tree(Spheres &spheres)
  : Sphere_Container( spheres ), _array( spheres ), _root( spheres.bad_sphere_index() ), _left(0), _right(0),
  _max_depth(0)
  { new_memory(); };

  virtual ~KD_Tree()
  { delete_memory(); }

  // Sphere_Container interface functions

  // get rid of contents of the current tree
  virtual void clear();

  // add one sphere to a tree, and the array
  virtual void add_sphere(size_t isphere);
  
  // after add_sphere, one may check whether the tree needs rebalancing
  bool needs_rebalance() const;

  // rebalance the tree for more uniform depth
  virtual void rebalance();

  // number of elements in the tree
  virtual size_t size() const {return _array.size();}

  // access the array form of the tree
  const Sphere_Container_Array &array() { return _array; }

  // tree[i]
  //
  // rvalue access only through operator []
  // size_t node = tree[9];
  const size_t & operator[](size_t i) const
  {
    assert(i < size());
    return _array[i];
  };
  // lvalue 
  size_t & operator[](std::size_t i)
  {
    return _array[i];      /* actual access, e.g. return mVector[idx]; */
  };

protected:
  // methods

  // add node to tree
  // the node has to be added to the array separately, e.g. by the search_tree_near_spheres method
  void add_sphere_to_tree(const size_t isphere);
  
  // memory allocation
  void new_memory();
  void delete_memory();

protected:

  // data

  // tree
  size_t _root;  
  size_t* _left; 
  size_t* _right; 

  // keep track of tree depth during insertions
  size_t _max_depth;

  // straight array of all the indices in the tree
  // takes space and time to build, but fast to access
  Sphere_Container_Array _array;
  
  virtual void reset_iterator_it( size_t &it ) const
  {
    return _array.reset_iterator_it( it );
  }
  virtual size_t next_it( size_t &it ) const
  {
    return _array.next_it( it );
  }

private:

  void clear_tree();

  // recursive functions called by public non-recursive function
  void clear_recurse(size_t i_sphere);
  void rebalance_recurse(size_t d, size_t lo, size_t hi);


};

#endif
