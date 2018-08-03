// Search_Grid.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0


#include "Search_Grid.hpp"


double Search_Grid::nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold )
{
  _stats.start_clock();

  Near_Workspace &ws = _ws;

  ws._pair_dist_threshold = ( distance_threshold < sqrt( std::numeric_limits<double>::max() ) ) ?
    distance_threshold * distance_threshold : std::numeric_limits<double>::max();

  // stuff that would normally be passed to the recursive function
  ws._conflict_sphere = Spheres::bad_sphere_index();
  
  point_to_cell_index( p, ws._current_cell);
  
  ws._p = p;

  ws._sub_container = 0;
  ws._num_calls = 0;

  nearest_sphere_recurse( 0 );

  near_sphere = ws._conflict_sphere;

  _stats.collect_stats(size(), ws._num_calls, 1);

  return sqrt(ws._pair_dist_threshold);
}


void Search_Grid::nearest_sphere_recurse( size_t d )
{
  Near_Workspace &ws = _ws;

  // debug
  assert(ws._num_calls < _num_cells * num_dim());
  ++ws._num_calls;
  
  const size_t next_d = d+1;
  if (next_d > num_dim())
  {
    // get nearest sphere in this cell
    Cell *cell = cell_from_cell_index(ws._current_cell);
    for (size_t j = 0; j < cell->size(); ++j)
    {
      const size_t sj = (*cell)[j];
      const double *s = _spheres[sj];
      if (ws._p != s)
      {
        const double pair_distance = _tt.distance_squared( ws._p, s);
        if (pair_distance <= ws._pair_dist_threshold)
        {
          ws._pair_dist_threshold = pair_distance;
          ws._conflict_sphere = sj;
        }
      }
    }
    return;
  }
  
  // recurse by dimension
  nearest_sphere_recurse( next_d );

  // iterate over the smaller and bigger cell indices, and recurse on them
  int original_index = ws._current_cell[d];
  int &new_index = ws._current_cell[d]; // alias
  double *cell_point = _tt.new_point();
  // for increment = -1 and +1
  for ( int increment = -1; increment < 2; increment +=2 )
  {
    new_index = original_index;
    while (1)
    {
      new_index += increment;
      
      // quit if moved from inside domain to outside
      if ((increment > 0 && new_index >= _axis_divisions) || (increment < 0 && new_index < 0))
        break;

      // is the next cell close enough to possibly contain the closest point?
      closest_point( ws._current_cell, ws._p, cell_point );
      const double pair_distance = _tt.distance_squared( ws._p, cell_point);
      // extending this next check to the power distance will depend on L,
      //  since we need some limit on the radius of a sphere in the next cell
      if (pair_distance <= ws._pair_dist_threshold)
      {
        // still outside the domain, don't recurse on the box but keep incrementing the index
        if (!((increment < 0 && new_index >= _axis_divisions) || (increment > 0 && new_index < 0)))
        {
          // recurse on neighboring cell
          nearest_sphere_recurse( next_d );
        }
      }
      else
        // cells are too far away to contain the closest point, subsequent cells will be even farther
        break;
    }
  }
  
  // clean_up:
  _tt.delete_point(cell_point);
  ws._current_cell[d] = original_index;
}



void Search_Grid::all_near_spheres(Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius )
{
  _stats.start_clock();
  sub_container->clear();

  // stuff that would normally be passed to the recursive function
  Near_Workspace &ws = _ws;
  
  point_to_cell_index( p, ws._current_cell);
  
  ws._p = p;
  double thresh = dist;
  if (subtract_radius)
    thresh += _spheres.radius(0);
  ws._pair_dist_threshold = thresh * thresh;

  ws._sub_container = sub_container;
  ws._num_calls = 0;
  ws._distances_checked = 0;

  all_near_spheres_recurse( 0 );

  
  _stats.collect_stats(size(), ws._distances_checked, sub_container->size());
}

void Search_Grid::all_near_spheres_recurse( size_t d )
{
  Near_Workspace &ws = _ws;

  // debug
  assert(ws._num_calls <= _num_cells * num_dim());
  ++ws._num_calls;
  
  const size_t next_d = d+1;
  if (next_d > num_dim())
  {
    // get near spheres in this cell
    Cell *cell = cell_from_cell_index(ws._current_cell);
    for (size_t j = 0; j < cell->size(); ++j)
    {
      const size_t sj = (*cell)[j];
      const double *s = _spheres[sj];
      if (ws._p != s)
      {
        const double pair_distance = _tt.distance_squared( ws._p, s);
        if (pair_distance < ws._pair_dist_threshold)
        {
          ws._sub_container->add_sphere(sj);
        }
      }
    }
    ws._distances_checked += cell->size();
    return;
  }
  
  // recurse by dimension
  int original_index = ws._current_cell[d];
  if ( original_index >= 0 && original_index < _axis_divisions )
    all_near_spheres_recurse( next_d );

  // iterate over the smaller and bigger cell indices, and recurse on them
  int &new_index = ws._current_cell[d]; // alias
  double *cell_point = _tt.new_point();
  // for increment = -1 and +1
  for ( int increment = -1; increment < 2; increment +=2 )
  {
    new_index = original_index;
    while (1)
    {
      new_index += increment;
      
      // quit if moved from inside domain to outside
      if ((increment > 0 && new_index >= _axis_divisions) || (increment < 0 && new_index < 0))
        break;

      // is the next cell close enough to possibly contain a near point?
      closest_point( ws._current_cell, ws._p, cell_point );
      const double pair_distance = _tt.distance_squared( ws._p, cell_point);
      // this fixed r isn't right. should depend on the cell contents, L or something
      if (pair_distance < ws._pair_dist_threshold)
      {
        // still outside the domain, don't recurse on the box but keep incrementing the index
        if (!((increment < 0 && new_index >= _axis_divisions) || (increment > 0 && new_index < 0)))
        {
          // recurse on neighboring cell
          all_near_spheres_recurse( next_d );
        }
      }
      else
        // cells are too far away to contain a close point, subsequent cells will be even farther
        break;
    }
  }
  
  // clean_up:
  _tt.delete_point(cell_point);
  ws._current_cell[d] = original_index;
}



bool Search_Grid::no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor, const double dist)
{
  _stats.start_clock();
  // stuff that would normally be passed to the recursive function
  Near_Workspace &ws = _ws;

  point_to_cell_index( p, ws._current_cell);
  
  ws._p = p;
  const double thresh = radius_factor * _spheres.radius(0)  + dist - 1e-5;
  ws._pair_dist_threshold = thresh*thresh;
  
  ws._sub_container = 0;
  ws._num_calls = 0;
  
  ws._conflict_sphere = _spheres.bad_sphere_index();
  
  bool return_value = no_near_spheres_recurse(0);
  near_sphere = ws._conflict_sphere;
  
  _stats.collect_stats(size(), ws._num_calls, !return_value);

  return return_value;
}
bool Search_Grid::no_near_spheres_recurse( size_t d )
{

  Near_Workspace &ws = _ws;

  // debug
  assert(ws._num_calls < (10+size()) * num_dim());
  ++ws._num_calls;

  // recurse by dimension
  const size_t next_d = d+1;
  // fathom
  // check spheres in the current cell
  if (next_d > num_dim())
  {
    Cell *cell = cell_from_cell_index(ws._current_cell);
    for (size_t j = 0; j < cell->size(); ++j)
    {
      const size_t sj = (*cell)[j];
      const double *s = _spheres[sj];
      if (ws._p != s)
      {
        const double pair_distance = _tt.distance_squared( ws._p, s);
        if (pair_distance < ws._pair_dist_threshold) // on surface should be OK
        {
          ws._conflict_sphere = sj;
          return false;
        }
      }
    }
    return true;
  }

  // recurse on current cell index
  if (!no_near_spheres_recurse( next_d ))
  {
    return false;
  }

  // iterate over the smaller and bigger cell indices
  int original_index = ws._current_cell[d];
  int &new_index = ws._current_cell[d]; // alias
  bool ret_val = true;
  double *cell_point = _tt.new_point();
  // for increment = -1 and +1
  for ( int increment = -1; increment < 2; increment +=2 )
  {
    new_index = original_index;
    while(1)
    {
      new_index += increment;
      if (new_index >= _axis_divisions || new_index <= 0)
        // end of domain
        break;
      // is the next cell close enough to possibly contain a near point?
      closest_point( ws._current_cell, ws._p, cell_point );
      const double pair_distance = _tt.distance_squared( ws._p, cell_point);
      // this fixed r isn't right. should depend on the cell contents, L or something
      if (pair_distance < ws._pair_dist_threshold)
      {
        // next cell contains a near point?
        if (!no_near_spheres_recurse( next_d ))
        {
          ret_val = false;
          break;
        }
      }
      else
        // cells are too far away to contain a close point
        break;
    }
  }
  _tt.delete_point(cell_point);
  ws._current_cell[d] = original_index;
  return ret_val;
}


void Search_Grid::closest_point( Cell_Index ci, const double *p, double *q )
{
  for (size_t d=0; d < num_dim(); ++d)
  {
    const double lc = lower_coordinate(ci, d);
    if ( lc > p[d])
      q[d] = lc;
    else
    {
      const double uc = upper_coordinate(ci, d);
      if ( uc < p[d])
        q[d] = uc;
      else
        q[d] = p[d];
    }
  }
}

void Search_Grid::trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2, const double radius_factor  )
{
  assert(0); // this is not implemented yet. The following simply iterates over *all* the spheres and ignores the grid.
  _stats.start_clock();
  size_t true_trims(0);
  for (size_t i = _array.first(); i != _array.bad_sphere_index(); i = _array.next())
  {
    if (_tt.anchored_trim( c, u, spheres()[i], A, A_1, A_2, radius_factor ) )
      ++true_trims;
    assert( A_1 < A_2 ); // but either could be negative
    assert( fabs(A_1) > r ); // the segment must be outside the initial sphere
    assert( fabs(A_2) > r );
  }
  _stats.collect_stats(sphere_container()->size(), sphere_container()->size(), true_trims);
}
