// Search_Tree.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Search_Tree.hpp"

// singleton class
// variables that remain constant during the recursions for near points
// put in the class for reduced memory allocation or stack size within the recursion
class Near_Workspace
{
public:
  // data
  const double * _p;
  double _thresh, _thresh_squared;
  Sphere_Container *_sub_container;
  size_t _conflict_sphere;
  size_t _num_calls; // debug, search stats

  // how to get the data  
  static Near_Workspace *instance() 
  { 
    if (!_the_instance)
      _the_instance = new Near_Workspace();
    return _the_instance;
  };

  ~Near_Workspace() {}

private:
  // private constructor to prevent anyone for making these, as they are expensive
  Near_Workspace() {}

  static Near_Workspace *_the_instance;
};

Near_Workspace *Near_Workspace::_the_instance(0);

double Search_Tree::nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold )
{
  _stats.start_clock();

  Near_Workspace &ws = *Near_Workspace::instance();
  ws._p = p;
  ws._thresh = distance_threshold;
  ws._thresh_squared = ( distance_threshold < sqrt( std::numeric_limits<double>::max() ) ) ? 
    distance_threshold * distance_threshold : std::numeric_limits<double>::max();
  ws._sub_container = 0;

  ws._sub_container = 0;
  ws._num_calls = 0;

  nearest_sphere_recurse( 0, _root );

  near_sphere = ws._conflict_sphere;

  _stats.collect_stats(size(), ws._num_calls, 1);

  return ws._thresh;
}
void Search_Tree::nearest_sphere_recurse( size_t d_index, size_t node )
{
  Near_Workspace &ws = *Near_Workspace::instance();

  // debug
  assert(ws._num_calls < size());
  ++ws._num_calls;
  
  const double * p = ws._p;

  if (_spheres[node] != p)
  {
    const double pair_distance = _tt.distance_squared( p, _spheres[node]);

    // new nearest? update the threshold distance for bounding branches
    if (pair_distance <= ws._thresh_squared)
    {
      ws._thresh_squared = pair_distance;
      // thresh, not thresh_squared, is used for the range search, need to take sqrt now
      ws._thresh = sqrt( pair_distance ); 
      ws._conflict_sphere = node;
    }  
  }

  // recurse
  const size_t next_d = (d_index + 1 == num_dim()) ? 0 : d_index + 1;
  // visit right and left if the current node is close enough that 
  // the right or left neighbors may be close enough.
  if ( (_right[node] != _spheres.bad_sphere_index()) && ( p[d_index] + ws._thresh > _spheres[node][d_index]))
    nearest_sphere_recurse( next_d, _right[node]);
  if ( ( _left[node] != _spheres.bad_sphere_index()) && ( p[d_index] - ws._thresh < _spheres[node][d_index]))
    nearest_sphere_recurse( next_d,  _left[node]);

}


void Search_Tree::all_near_spheres(Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius )
{
  _stats.start_clock();
  sub_container->clear();
  if (_root == _spheres.bad_sphere_index())
    return;

  // stuff that would normally be passed to the recursive function
  Near_Workspace &ws = *Near_Workspace::instance();
  ws._p = p;
  ws._thresh = dist;
  if (subtract_radius)
    ws._thresh += _spheres.radius(0);
  ws._thresh_squared = ws._thresh * ws._thresh;

  ws._sub_container = sub_container;
  ws._num_calls = 0;

  
  // this function strongly affects the overall run-time
  all_near_spheres_recurse( 0, _root );

  _stats.collect_stats(size(), ws._num_calls, sub_container->size());
}
void Search_Tree::all_near_spheres_recurse( size_t d_index, size_t node )
{
  Near_Workspace &ws = *Near_Workspace::instance();

  // debug
  assert(ws._num_calls < size());
  ++ws._num_calls;
  
  const double * p = ws._p;

  if (_spheres[node] != p)
  {
    const double pair_distance = _tt.distance_squared( p, _spheres[node]);

    // close, so add to array
    if (pair_distance < ws._thresh_squared)
    {
      ws._sub_container->add_sphere(node);
    }  
  }

  // recurse
  const size_t next_d = (d_index + 1 == num_dim()) ? 0 : d_index + 1;
  // visit right and left if the current node is close enough that 
  // the right or left neighbors may be close enough.
  if ( (_right[node] != _spheres.bad_sphere_index()) && ( p[d_index] + ws._thresh > _spheres[node][d_index]))
    all_near_spheres_recurse( next_d, _right[node]);
  if ( ( _left[node] != _spheres.bad_sphere_index()) && ( p[d_index] - ws._thresh < _spheres[node][d_index]))
    all_near_spheres_recurse( next_d,  _left[node]);

}



bool Search_Tree::no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor, const double dist)
{
  _stats.start_clock();
  if (_root == _spheres.bad_sphere_index())
  {
    near_sphere = _spheres.bad_sphere_index();
    return true;
  }

  // stuff that would normally be passed to the recursive function
  Near_Workspace &ws = *Near_Workspace::instance();
  ws._p = p;
  ws._thresh = radius_factor * _spheres.radius(0)  + dist - 1e-5;
  ws._thresh_squared = ws._thresh * ws._thresh;

  ws._sub_container = 0;
  ws._num_calls = 0;

  ws._conflict_sphere = _spheres.bad_sphere_index();

  bool return_value = no_near_spheres_recurse( 0, _root );
  near_sphere = ws._conflict_sphere;
  
  _stats.collect_stats(size(), ws._num_calls, !return_value);

  return return_value;
}
bool Search_Tree::no_near_spheres_recurse( size_t d_index, size_t node ) 
{

  Near_Workspace &ws = *Near_Workspace::instance();

  // debug
  assert(ws._num_calls < size());
  ++ws._num_calls;

  if (ws._p != _spheres[node])
  {
    const double pair_distance = _tt.distance_squared( ws._p, _spheres[node]);

    if (pair_distance < ws._thresh_squared) // on surface should be OK
    {
      ws._conflict_sphere = node;
      return false; 
    }
  }

  // recurse
  const bool check_right = ((_right[node] != _spheres.bad_sphere_index()) && (ws._p[d_index] + ws._thresh >  _spheres[node][d_index]));
  const bool check_left  = (( _left[node] != _spheres.bad_sphere_index()) && (ws._p[d_index] - ws._thresh <  _spheres[node][d_index]));
  const size_t next_d = (d_index + 1 == num_dim()) ? 0 : d_index + 1;
  if (check_right && !no_near_spheres_recurse(next_d, _right[node])) return false;
  if (check_left  && !no_near_spheres_recurse(next_d,  _left[node])) return false;

  return true;
}

void Search_Tree::trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2, const double radius_factor )
{
  _stats.start_clock();
  if (_root == _spheres.bad_sphere_index())
    return;
  Trim_Statespace * z = &trim_statespace;
  z->set(c, u, r, A, A_1, A_2, radius_factor);
  z->_visited = 0; z->_found = 0;
  trim_line_anchored_recurse( 0, _root );
  _stats.collect_stats(size(), z->_visited, z->_found);
}
void Search_Tree::trim_line_anchored_recurse( size_t d_index, size_t sphere_index )                                                    
{
  // indices to traverse the kd-tree
  // size_t d_index, size_t sphere_index
  assert(sphere_index != _spheres.bad_sphere_index());

  if (d_index == num_dim()) d_index = 0;
  const double *s = _spheres[sphere_index];

  // the geometry work
  Trim_Statespace * z = &trim_statespace;
  z->anchored_trim(_tt, s);
  
  const double &x_s = s[d_index]; // "x" coord of test sphere

  // this right and left test assumes that z.r is the radius of any future s
  // could be modified to be a maximum...

  // right?
  {
    const double x_1 = z->_c[d_index] + *(z->_A_1) * z->_u[d_index];   // "x" coord of point at A_1
    const double x_2 = z->_c[d_index] + *(z->_A_2) * z->_u[d_index];   // "x" coord of point at A_2
    const size_t &next_s = _right[sphere_index];
    if ((next_s != _spheres.bad_sphere_index()) && ( x_1 > x_s - z->_r || x_2 > x_s - z->_r ) )
    {
      trim_line_anchored_recurse(d_index + 1, next_s);
    }
  }

  // left? 
  // re-do x_1 and x_2, since the right recursion might have trimmed these
  {
    const double x_1 = z->_c[d_index] + *(z->_A_1) * z->_u[d_index];
    const double x_2 = z->_c[d_index] + *(z->_A_2) * z->_u[d_index];
    const size_t &next_s = _left[sphere_index];
    if ((next_s != _spheres.bad_sphere_index()) && ( x_1 < x_s + z->_r || x_2 < x_s + z->_r ))
    {
      trim_line_anchored_recurse(d_index + 1, next_s);
    }
  }

}
