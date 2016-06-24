// Search_Ghost_Structure.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef SEARCH_GHOST_STRUCTURE_HPP
#define SEARCH_GHOST_STRUCTURE_HPP

// generic interface for searching over ghosts

#include "Search_Structure.hpp"
#include "Sphere_Container_Array.hpp"
#include "Ghost_Spheres.hpp"
#include "Ghost_Global_Container.hpp"

class Search_Ghost_Structure : public Search_Structure
{
public:
  // ================================================
  // Search_Ghost specific
  // ================================================
  Search_Ghost_Structure( Ghost_Spheres &ghost_spheres )
  : Search_Structure( ghost_spheres.num_dim() ),
  _ghost_spheres(ghost_spheres), _search(0),
  _search_results( ghost_spheres ), _all_spheres(ghost_spheres), _reghost(true)
  {}

  Search_Ghost_Structure( Ghost_Spheres &ghost_spheres, Search_Structure *search )
  : Search_Structure( ghost_spheres.num_dim() ),
  _ghost_spheres(ghost_spheres), _search(search),
  _search_results( ghost_spheres ), _all_spheres(ghost_spheres), _reghost(true)
  {}

  virtual ~Search_Ghost_Structure() {}

  // by default, the searches automatically instantiate the needed ghosts
  void set_search_reghost( bool reghost ) { _reghost = reghost; }

  // This class just adds a layer to handle ghosting; underlying it is a real search structure.
  // This might be a Search_Tree, or a Search_Array,
  // or might be implemented some other way.
  void set_underlying_search( Search_Structure *search ) {_search = search;}


  // ================================================
  // Overloaded interface from Search_Structure
  // ================================================

  // Warning: all_near_spheres may automatically dematerialize old ghosts and materializes new ones
  //   in the local neighborhood of the point. Also makes periodic copies.
  // For efficiency, subsequent searches should use the sub_container
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() ) = 0;

  virtual void all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false ) = 0;

  // Warning: does not make temporary ghosts, does make periodic copies of p.
  // the near_sphere may be near a periodic copy of p, rather than p itself
  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.) = 0;

  // not implemented yet. Would be easier if we had "all_near_spheres" for line segments
  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
    const double radius_factor = 1. )
  {
    assert( false ); // not implemented
    assert(_search);
    _search->trim_line_anchored(c, u, r, A, A_1, A_2, radius_factor );
  }

  virtual Sphere_Container *sphere_container() {return &_all_spheres;}
  virtual bool balance_container() { assert(_search); return _search->balance_container(); }
  virtual Rebalance_Stats *rebalance_stats() { assert(_search); return _search->rebalance_stats();}
  
  virtual void add_sphere( size_t si )
  { sphere_container()->add_sphere(si); _search->add_sphere(si); }


protected:
  Ghost_Spheres &_ghost_spheres;
  Search_Structure *_search;
  Sphere_Container_Array _search_results;
  bool _reghost;
  Ghost_Global_Container _all_spheres;
  
  // shorthand
  size_t first_real() { return advance_to_real( _search->sphere_container()->first() ); }
  size_t bad_sphere_index() { return Spheres::bad_sphere_index(); }
  size_t next_real() { return advance_to_real( _search->sphere_container()->next() ); }
  size_t advance_to_real( size_t i )
  {
    while ( (i != Spheres::bad_sphere_index()) && _ghost_spheres.is_ghost(i) )
    {
      i = _search->sphere_container()->next();
    }
    return i;
  }
  const size_t num_dim() { return _ghost_spheres.num_dim(); }


private:

};

#endif