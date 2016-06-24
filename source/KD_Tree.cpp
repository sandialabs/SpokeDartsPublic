// KD_Tree.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "KD_Tree.hpp"
#include <assert.h>

void KD_Tree::new_memory()
{
  // delete old memory, if any
  delete_memory();

  // allocate new
  size_t max_size = _spheres.max_num_spheres();
  if (max_size)
  {
    _left  = new size_t[max_size];
    _right = new size_t[max_size];
    for (size_t i = 0; i < max_size; i++)
    {
      _left[i] = _spheres.bad_sphere_index();
      _right[i] = _spheres.bad_sphere_index();
    }
  }
}

void KD_Tree::delete_memory()
{
  delete []  _left;  _left = 0;
  delete [] _right; _right = 0;
}

void KD_Tree::clear()

{
  clear_tree();

  _array.clear();
}
void KD_Tree::clear_tree()
{
  // traverse the tree and reset the pointers
  if ((_root != _spheres.bad_sphere_index()) && _left)
    clear_recurse(_root);
  _root = _spheres.bad_sphere_index();
  _max_depth = 0;
  
  // expensive debug test
//  if (false)
//  {
//    for (size_t i = 0; i < _spheres.size(); ++i)
//    {
//      assert( _left == 0 ||  _left[i] == _spheres.bad_sphere_index());
//      assert(_right == 0 || _right[i] == _spheres.bad_sphere_index());
//    }
//  }
}
void KD_Tree::clear_recurse(size_t isphere)
{
  // isphere might be a temporary ghost index that has been dematerialized already. That is OK here.
  assert( _spheres.is_valid_index( isphere ) );
  assert(_left);
  assert(_right);
  
  if (_left[isphere] != _spheres.bad_sphere_index())
  {
    clear_recurse(_left[isphere]);
    _left[isphere] = _spheres.bad_sphere_index();
  }
  if (_right[isphere] != _spheres.bad_sphere_index())
  {
    clear_recurse(_right[isphere]);
    _right[isphere] = _spheres.bad_sphere_index();
  }
}


void KD_Tree::rebalance()
{
  if ( !size() )
    return;
  
  // debug
  const double past_time = _rebalance_stats.cpu_time;
  
  _rebalance_stats.start_clock();
  _rebalance_stats.num_rebalances++;
  
  clear_tree(); // leave array intact, it is used below
  rebalance_recurse( 0, 0, size()-1 );

  _rebalance_stats.collect_stats();
  const double elapsed_time = _rebalance_stats.cpu_time - past_time;
  
  // std::cout << "Rebalance " << _rebalance_stats.num_rebalances << " took " << elapsed_time << " seconds." << std::endl;
}

struct xd_less
{
  xd_less(size_t d, Spheres &spheres) : _d(d), _spheres(spheres) {}
  size_t _d;
  Spheres &_spheres;
  inline bool operator() (const size_t &pi, const size_t &qj)
  {
    return ( _spheres[pi][_d] <  _spheres[qj][_d] );
  }
};


void KD_Tree::rebalance_recurse(size_t d, size_t lo, size_t hi)
{
  if (d > num_dim())
    d = 0;

  // sort on the d_th coordinate
  xd_less my_less(d, _spheres);
  size_t *lop = &_array[lo];
  size_t *hip = &_array[hi+1];
  std::sort( lop, hip, my_less );

  // find median. 
  // bal alternates putting the larger half in the right or the left
  const size_t bal = d % 2;
  size_t mid = (lo + hi + bal) / 2;

  add_sphere_to_tree( _array[mid] );

  if (mid < hi)
    rebalance_recurse( d+1, mid+1, hi );
  if (lo < mid )
    rebalance_recurse( d+1, lo, mid-1 );
}

void KD_Tree::add_sphere(const size_t isphere)
{  
  // add to array
  _array.add_sphere(isphere);

  add_sphere_to_tree( isphere );
}

void KD_Tree::add_sphere_to_tree(const size_t isphere)
{
  // add to tree

  // first node is the root
  if (_root == _spheres.bad_sphere_index())
  {
    _root = isphere;
    return;
  }

  assert( _spheres.is_valid_sphere( isphere ) );

  // insert sphere into tree
  size_t node(_root), d_index(0), depth(0);
  while (true)
  {
    depth++;
    const double & x_tree_sphere = _spheres[node][d_index];
    const double & x_sphere = _spheres[isphere][d_index];
    if (x_sphere >  x_tree_sphere)
    {
      if (_right[node] != _spheres.bad_sphere_index())
        node = _right[node];
      else
      {
        _right[node] = isphere;
        if (depth > _max_depth)
          _max_depth = depth;
        return;
      }
    }
    else
    {
      if (_left[node] != _spheres.bad_sphere_index())
        node = _left[node];
      else
      {
        _left[node] = isphere;
        if (depth > _max_depth)
          _max_depth = depth;
        return;
      }
    }

    if (++d_index == num_dim())
      d_index = 0;
  }
  // unreachable
}

bool KD_Tree::needs_rebalance() const
{
  if (size() < 16)
    return false;

  // how far out of balance are we?
  // perfect is log2( size )
  // ok is log size * log log size
  size_t sz = size();
  size_t mylog(0);
  while ( sz >>= 1 ) ++mylog;
  size_t m( mylog );
  size_t myloglog(0);
  while ( m >>= 1 ) ++myloglog;
  
  const size_t perfect_depth = mylog;
  const size_t ok_depth = perfect_depth * (size_t) ceil( sqrt(perfect_depth) ); // or * myloglog;
  // to do: also put a cap on how often rebalancing is done in the calling loop
  // debug
//  if (_max_depth > ok_depth)
//    std::cout << "Rebalancing depth:" << _max_depth << " size:" << size() << " perfect_depth: " << perfect_depth << " OK_depth:" << ok_depth << std::endl;
  return (_max_depth > ok_depth);
}