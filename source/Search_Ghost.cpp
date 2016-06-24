// Search_Ghost.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0


#include "Search_Ghost.hpp"


double Search_Ghost::nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold )
{
  // uses workspace data members _g and _offset
  _stats.start_clock();
  
  if (_reghost)
    _ghost_spheres.dematerialize_temporary_ghosts();
  const size_t num_searched = _ghost_spheres.size();
  
  // in every case, add the results of the original point, without ghosting
  // do this first, before any new local ghosts are materialized by the periodic copies
  assert(_search);
  const size_t start_visited = _search->_stats.visited();
  // update distance_threshold "dt" with the nearest found sphere so far
  double dt = _search->nearest_sphere( p, near_sphere, distance_threshold );
  
  // find all the real and materialized points
  // materialize more if needed
  // speed up, only do this if the point is within the negative frame of the domain.
  if ( _ghost_spheres.is_temporary_ghosts() && !_ghost_spheres.no_ghosts_needed(p, dt) )
  {
    // make periodic copies of the query point, and search in the vicinity of each
    // convert the nearby points to ghosts and return their indices in the sub_structure
    if (_reghost)
      _ghost_spheres.make_periodic_sphere_copies(p);
    const Sphere_Array& pc( _ghost_spheres.periodic_copies() ); // shorthand
    const size_t& pc_size( _ghost_spheres.periodic_copies_size() ); // shorthand
    _tt.radius(_g) = _tt.radius(p);
    bool orig_is_closest = true;
    size_t closest_i ( near_sphere );
    for (size_t i = 0; i < pc_size; ++i)
    {
      if ( !_ghost_spheres.is_orig_index(i) )
      {
        // find the nearest disk to the periodic copy of p
        size_t near_sphere_i;
        double dti = _search->nearest_sphere( pc[i], near_sphere_i, dt );
        if ( dti < dt )
        {
          if ( near_sphere_i != Spheres::bad_sphere_index() )
          {
            dt = dti;
            closest_i = i;
            orig_is_closest = false;
            near_sphere = near_sphere_i;
          }
        }
      }
    }
    
    // the nearest point is offset by the period of the domain
    // translate it into the frame of the original point
    // and update the near_sphere to point to this ghost
    if (!orig_is_closest)
    {
      _tt.subtract(_offset, p, pc[ closest_i ] );
      _tt.add(_g, _offset, _ghost_spheres[near_sphere] );
      near_sphere = _ghost_spheres.add_ghost_sphere_to_array(_g);
    }
  }
  
  const size_t visited = _search->_stats.visited() - start_visited;
  _stats.collect_stats( num_searched, visited, 1 );

  return dt;
}


void Search_Ghost::all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius )
{
  // uses workspace data members _g and _offset
  _stats.start_clock();
  
  if (_reghost)
    _ghost_spheres.dematerialize_temporary_ghosts();
  const size_t num_searched = _ghost_spheres.size();

  // in every case, add the results of the original point, without ghosting
  // do this first, before any new local ghosts are materialized by the periodic copies
  assert(_search);
  const size_t start_visited = _search->_stats.visited();
  _search->all_near_spheres( sub_container, p, dist, subtract_radius );

  // find all the real and materialized points
  // materialize more if needed
  // speed up, only do this if the point is within the negative frame of the domain.
  if ( _ghost_spheres.is_temporary_ghosts() && !_ghost_spheres.no_ghosts_needed(p, dist) )
  {
    // make periodic copies of the query point, and search in the vicinity of each
    // convert the nearby points to ghosts and return their indices in the sub_structure
    if (_reghost)
      _ghost_spheres.make_periodic_sphere_copies(p);
    const Sphere_Array& pc( _ghost_spheres.periodic_copies() ); // shorthand
    const size_t& pc_size( _ghost_spheres.periodic_copies_size() ); // shorthand
    _tt.radius(_g) = _tt.radius(p);
    for (size_t i = 0; i < pc_size; ++i)
    {
      if ( !_ghost_spheres.is_orig_index(i) )
      {
        // all the disks near the periodic copy of p
        _search->all_near_spheres( &_search_results, pc[i], dist, subtract_radius );
        
        // the found points are offset by the period of the domain
        // translate them into the frame of the original point
        _tt.subtract(_offset, p, pc[i] );

        // add the found-and-offset points as ghosts and to the subcontainer
        for (size_t ii = _search_results.first(); ii != _ghost_spheres.bad_sphere_index(); ii = _search_results.next() )
        {
          assert( ii != Spheres::bad_sphere_index() );
          assert( ii < _ghost_spheres.Spheres::size() ); // i.e. not a ghost
          _tt.add(_g, _offset, _ghost_spheres[ii] );
          // create a temporary ghost
          size_t gi = _ghost_spheres.add_ghost_sphere_to_array(_g);
          // add the ghost to the sub_container
          sub_container->add_sphere(gi);
          // double-check that g's distance to p is less than dist
          assert( _tt.distance_squared(_g,p) - (subtract_radius ? _tt.radius(_g)*_tt.radius(_g) : 0.)  < (dist * dist)*1.01 );
        }
      }
    }
  }

  const size_t found = sub_container->size();
  const size_t visited = _search->_stats.visited() - start_visited;
  assert( visited >= found );
  _stats.collect_stats( num_searched, visited, found );
}

bool Search_Ghost::no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor, const double dist)
{
  _stats.start_clock();

  const size_t start_visited = _search->_stats.visited();

  // check original nbhd
  if (!_search->no_near_spheres( p, near_sphere, radius_factor, dist) )
  {
    const size_t visited = _search->_stats.visited() - start_visited;
    _stats.collect_stats( _ghost_spheres.size(), visited, 1 );
    return false;
  }

  // find all the real and materialized points
  if ( _ghost_spheres.is_temporary_ghosts())
  {
    // make periodic copies of the query point, and search in the vicinity of each
    if (_reghost)
      _ghost_spheres.make_periodic_point_copies(p);
    const Sphere_Array& pc( _ghost_spheres.periodic_copies() ); // shorthand
    const size_t& pc_size( _ghost_spheres.periodic_copies_size() ); // shorthand
    for (size_t i = 0; i < pc_size; ++i)
    {
      if ( !_ghost_spheres.is_orig_index(i) )
      {
        if (!_search->no_near_spheres( pc[i], near_sphere, radius_factor, dist) )
        {
          const size_t visited = _search->_stats.visited() - start_visited;
          _stats.collect_stats( _ghost_spheres.size(), visited, 1 );
          return false;
        }
      }
    }
  }
  const size_t visited = _search->_stats.visited() - start_visited;
  _stats.collect_stats( _ghost_spheres.size(), visited, 0 );
  return true;
}
