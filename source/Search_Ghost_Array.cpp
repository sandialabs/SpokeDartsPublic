// Search_Ghost_Array.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0


#include "Search_Ghost_Array.hpp"

// find the ghost of q that is closest to p, put it in _g


double Search_Ghost_Array::nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold )
{
  // uses workspace data members _g and _offset
  _stats.start_clock();
  
  if (_reghost)
    _ghost_spheres.dematerialize_temporary_ghosts();
  const size_t num_searched = _ghost_spheres.size();

  assert(_search);

  double thresh_squared = ( distance_threshold < sqrt( std::numeric_limits<double>::max() ) ) ?
  distance_threshold * distance_threshold : std::numeric_limits<double>::max();
  near_sphere = bad_sphere_index();
  
  // go over all the real spheres
  for (size_t i = first_real(); i != bad_sphere_index(); i = next_real())
  {
    const double *q = _ghost_spheres[i];
    if (q != p)
    {
      // find the ghost that is closest, sets _g
      double dx2 = _tt.closest_ghost_distance_squared(p, q);

      if ( dx2 < thresh_squared)
      {
        thresh_squared = dx2;
        near_sphere = i;
      }
    }
  }

  // none close enough found
  if ( near_sphere == bad_sphere_index() )
    return distance_threshold;
  
  //
  bool is_q(true);
  _tt.closest_ghost(p, _ghost_spheres[near_sphere], _g, is_q);
  if ( !is_q )
  {
    near_sphere = _ghost_spheres.add_ghost_sphere_to_array( _g );
  }
  _stats.collect_stats(num_searched, num_searched, 1);
  return sqrt(thresh_squared);
}


void Search_Ghost_Array::all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius )
{
  // uses workspace data members _g and _best_g
  _stats.start_clock();
  
  sub_container->clear();
  
  if (_reghost)
    _ghost_spheres.dematerialize_temporary_ghosts();
  const size_t num_searched = _ghost_spheres.size();

  assert(_search);
  
  const double dist_squared = ( dist < sqrt( std::numeric_limits<double>::max() ) ) ?
  dist * dist : std::numeric_limits<double>::max();
  
  // go over all the real spheres
  bool is_q(true);
  for (size_t i = first_real(); i != bad_sphere_index(); i = next_real())
  {
    const double *q = _ghost_spheres[i];
    if (q != p)
    {
      double dx2 = _tt.closest_ghost(p, q, _g, is_q);
      
      if (subtract_radius)
        dx2 -= _tt.radius(q) * _tt.radius(q);
        
      if ( dx2 < dist_squared)
      {
        sub_container->add_sphere( is_q ? i:
                                  _ghost_spheres.add_ghost_sphere_to_array( _g ));
      }
    }
  }

  const size_t found = sub_container->size();
  _stats.collect_stats( num_searched, num_searched, found );
}

bool Search_Ghost_Array::no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor, const double dist)
{
  _stats.start_clock();

  if (_reghost)
    _ghost_spheres.dematerialize_temporary_ghosts();
  const size_t num_searched = _ghost_spheres.size();

  // tolerance will skip some almost-close spheres
  const double thresh = dist + _ghost_spheres.radius(0) * radius_factor - 1E-10;
  const double thresh_squared = thresh * thresh;
  
  // go over all the real spheres, quitting if I find a close one or a close ghost
  for (size_t i = first_real(); i != bad_sphere_index(); i = next_real())
  {
    const double *q = _ghost_spheres[i];
    if (q != p)
    {
      // find the ghost that is closest, sets _g
      double dx2 = _tt.closest_ghost_distance_squared(p, q);
      
      if ( dx2 < thresh_squared)
      {
        _stats.collect_stats(num_searched, i, 1);
        return false;
      }
    }
  }

  _stats.collect_stats( num_searched, num_searched, 0 );
  return true;
}

void Search_Ghost_Array::trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
                        const double radius_factor )
{
  _stats.start_clock();
  size_t true_trims(0);
  
  double a = A;
  const double max_dist = 2. * num_dim(); // big number, more than box diagonal
  double second_dist2( max_dist ), min_second_dist2( max_dist );
  
  // first pass: go over all the points in sequence
  
  // convert anchor point to coordinates
  _tt.axpy(_anchor, a, u, c );
  bool is_q(false);
    
  // trim by the closest ghost to the anchor
  for (size_t i = first_real(); i != bad_sphere_index(); i = next_real())
  {
    // set _g to be the closest ghost to the marching point a
    const double *q = _ghost_spheres[i];
    _tt.two_closest_ghosts(_anchor, q, _g, is_q, second_dist2);
    if ( second_dist2 < min_second_dist2 )
      min_second_dist2 = second_dist2;
    
    if (_tt.anchored_trim( c, u, _g, A, A_1, A_2, radius_factor ) )
      ++true_trims;
    assert( A_1 <= A_2 ); // but either could be negative
  }
  
  // second pass: march each point using its ghosts to A_2
  const double dA = sqrt(min_second_dist2) - r; // radius_factor already included in r
  
  // for extremely large radius, the second closest ghost might already cover the anchor!
  // but that is suposed to be caught before this function is called
  assert( dA > 0 );
  
  a = A + dA;
  if ( a < A_2 )
  {
    bool giveup(false);
    for (size_t i = first_real(); i != bad_sphere_index() && !giveup; i = next_real())
    {
      const double *q = _ghost_spheres[i];
      double aa = a;
      if ( aa >= A_2 )
        break; // quit looping over the points, as we've trimmed it far back already
      // march the anchor until A_2 is passed
      // note A_2 can change
      size_t itercount = 0;
      const size_t maxiter = 100;
      size_t giveup_count(0);
      const size_t giveup_max(3);
      do {
        ++itercount;
        
        // march the anchor
        _tt.axpy(_anchor, aa, u, c );

        // set _g to be the closest ghost to the marching anchor
        _tt.two_closest_ghosts(_anchor, q, _g, is_q, second_dist2);
      
        // trim by g
        // use the true A here
        if (_tt.anchored_trim( c, u, _g, A, A_1, A_2, radius_factor ) )
          ++true_trims;
        assert( A_1 <= A_2 ); // but either could be negative
        
        // march aa
        const double daa = sqrt( second_dist2 ) - r; // radius_factor included already
        assert( daa > 0. );
        
        // give up if stepsize is tiny or too many iterations
        if ( (fabs(daa) < 1.e-3 * _tt.radius(q)) || (itercount > maxiter) )
        {
          ++giveup_count;
          giveup = giveup_count > giveup_max;
        }
        // we can only count on the spoke being uncovered up to aa
        if (giveup && (A_2 > aa))
          A_2 = aa; // max we can trust, before marching
        
        aa += daa;
      } while ( aa < A_2 && !giveup);
    }
  }
  
  // third pass: march each point using its ghosts to A_1, if needed
  a = A - dA;
  if ( a > A_1 )
  {
    bool giveup(false);
    for (size_t i = first_real(); i != bad_sphere_index(); i = next_real())
    {
      const double *q = _ghost_spheres[i];
      double aa = a;
      if ( aa <= A_1 )
        break; // quit looping over the points, as we've trimmed it far back already
      // march the anchor until A_1 is passed
      // note A_1 can change
      size_t itercount = 0;
      const size_t maxiter = 100;
      size_t giveup_count(0);
      const size_t giveup_max(3);
      do {
        ++itercount;

        // march the anchor
        _tt.axpy(_anchor, aa, u, c );
        
        // set _g to be the closest ghost to the marching anchor
        _tt.two_closest_ghosts(_anchor, q, _g, is_q, second_dist2);
        
        // trim by g
        // use the true A here
        if (_tt.anchored_trim( c, u, _g, A, A_1, A_2, radius_factor ) )
          ++true_trims;
        assert( A_1 <= A_2 ); // but either could be negative
        
        const double daa = sqrt( second_dist2 ) - r; // radius_factor included already
        assert( daa > 0. );

        // give up if stepsize is tiny or too many iterations
        if ( (fabs(daa) < 1.e-3 * _tt.radius(q)) || (itercount > maxiter) )
        {
          ++giveup_count;
          giveup = giveup_count > giveup_max;
        }
        // we can only count on the spoke being uncovered up to aa
        if (giveup && (A_1 < aa))
          A_1 = aa; // max we can trust, before marching

        // march aa
        aa -= daa;
      
      } while ( aa > A_1 && !giveup );
    }
  }  
  _stats.collect_stats(sphere_container()->size(), sphere_container()->size(), true_trims);
}

