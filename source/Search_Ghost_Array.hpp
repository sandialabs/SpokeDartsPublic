// Search_Ghost_Array.hpp


// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0


// method that avoids exponential complexity in dimension

#ifndef SEARCH_GHOST_ARRAY_HPP
#define SEARCH_GHOST_ARRAY_HPP

#include "Search_Ghost_Structure.hpp"
#include "Search_Array.hpp"

class Search_Ghost_Array : public Search_Ghost_Structure
{
public:
  // ================================================
  // Search_Ghost_Array specific
  // ================================================
  Search_Ghost_Array( Ghost_Spheres &ghost_spheres, Search_Array *search )
  : Search_Ghost_Structure( ghost_spheres, search ),
  _best_g( _tt.new_sphere() ), _g( _tt.new_sphere() ), _anchor( _tt.new_sphere() )
  {}

  ~Search_Ghost_Array()
  { _tt.delete_sphere(_best_g); _tt.delete_sphere(_g); _tt.delete_sphere(_anchor); }

  // ================================================
  // Overloaded interface from Search_Structure
  // ================================================

  // Warning: all_near_spheres may automatically dematerialize old ghosts and materializes new ones
  //   in the local neighborhood of the point. Also makes periodic copies.
  // For efficiency, subsequent searches should use the sub_container
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() );

  virtual void all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false );

  // The near_sphere may be near a periodic copy of p, rather than p itself
  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.);

  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
                                  const double radius_factor = 1. );
  

protected:
  
private:
  // workspace for all_near_spheres
  double *_g, *_best_g, *_anchor;
  

};

#endif