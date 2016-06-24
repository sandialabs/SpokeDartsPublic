// Search_Ghost.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef SEARCH_GHOST_HPP
#define SEARCH_GHOST_HPP

#include "Search_Ghost_Structure.hpp"

class Search_Ghost : public Search_Ghost_Structure
{
public:
  // ================================================
  // Search_Ghost specific
  // ================================================
  Search_Ghost( Ghost_Spheres &ghost_spheres )
  : Search_Ghost_Structure( ghost_spheres ),
  _offset( _tt.new_sphere() ), _g( _tt.new_sphere() )
  {}

  Search_Ghost( Ghost_Spheres &ghost_spheres, Search_Structure *search )
  : Search_Ghost_Structure( ghost_spheres, search ),
  _offset( _tt.new_sphere() ), _g( _tt.new_sphere() )
  {}

  ~Search_Ghost() { _tt.delete_sphere(_offset); _tt.delete_sphere(_g); }


  // ================================================
  // Overloaded interface from Search_Structure
  // ================================================

  // Warning: all_near_spheres may automatically dematerialize old ghosts and materializes new ones
  //   in the local neighborhood of the point. Also makes periodic copies.
  // For efficiency, subsequent searches should use the sub_container
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() );

  virtual void all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false );

  // Warning: does not make temporary ghosts, does make periodic copies of p.
  // the near_sphere may be near a periodic copy of p, rather than p itself
  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.);

  // not implemented yet. Would be easier if we had "all_near_spheres" for line segments
  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
    const double radius_factor = 1. )
  {
    assert( false ); // not implemented
    assert(_search);
    _search->trim_line_anchored(c, u, r, A, A_1, A_2, radius_factor );
  }


protected:

private:
  // workspace for all_near_spheres
  double *_g, *_offset;


};

#endif