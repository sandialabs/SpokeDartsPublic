// Search_Structure.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// This defines the interface for standard searches, be they implemented by trees or arrays
// This is an abstract class; derived classes must fill in the implementations

#ifndef SEARCH_STRUCTURE
#define SEARCH_STRUCTURE

#include <cstddef>
#include <limits>
#include "Sphere_Container.hpp"
#include "Search_Stats.hpp"
#include "Trimming_Tool.hpp"

class Search_Structure
{
public:
  bool is_global() {return _is_global;}
  
  // return the distance to the nearest sphere to point p, excluding p itself
  // todo future: provide a variation that uses the power distance to the sphere
  // Only include spheres closer or equal to the distance_threshold: 
  //   return double-max and bad_sphere_index if there is no such nearest sphere.
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() ) = 0;

  // return true if no node (center + radius * radius_factor) is within distace to p
  //   does not check p itself
  // replaces kd_tree_valid_dart_neighbors_tree
  // replaces kd_tree_duplicate_intersection_point
  // replaces kd_tree_covered_intersection_points
  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.) = 0;

  // return, in the sub_container, all spheres near to the query spheres, and only those.
  //   does not include p itself
  // (The container is emptied first.)
  // near = spheres s (center) are within distance dist or dist-radius(s) from the point p
  // was kd_tree_retrieve_neighbors
  virtual void all_near_spheres(  Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false ) = 0;


  // in the future, we may want a line-segment version of no_near_spheres and all_near_spheres

  // if you want an un-anchored trim, then use piercings
  // it doesn't use the tree to gather, and just intersects all the spheres passed in
  // but you can use "all_near_spheres" for segments to find a subtree to search instead

  // anchored trim, just retain the one uncovered line segment containing the anchor
  // replaces kd_tree_trim_anchored_line_spoke
  // replaces kd_tree_trim_line_spoke
  // currently, for the tree traversal to work, all the spheres are of radius less than r,
  // and r must be the radius_factor * the radius of any sphere
  // For a version that considers all segments, not just the anchored ones, see SpokeDarts::pick_segment_piercing
  // in the future, we may want an axis-aligned version, for k-d darts
  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
    const double radius_factor = 1. ) = 0;
    // points or vectors
    // p -> c  // center of sphere
    // q ->    // dart point
    //      u  // unit vector from c towards a, the spoke direcion
    //
    // distances, non-negative
    // r -> r  sphere radius
    //          A  distance from c to a
    //          A_1  spoke extends from A_1 to A along unit vector
    // t_end -> A_2  spoke extends from A to A_2
    //      : A_1 < A < A_2
    //
    // indices to traverse the kd-tree
    // size_t d_index, size_t sphere_index
    // if !use_extext_neighbors then use_ext_neighbors for the trimming

  // find uncovered arcs of a great circle
  // replaces gather_crossings, plus caller
  // For a version that considers all segments, not just the anchored ones, see SpokeDarts::wheel and SpokeDarts::pick_segment_crossing
  // virtual void trim_circle_anchored(); //zzyk to do

  // access the sphere_container that the search is performed over.
  // e.g. it may be needed for a call to all_near_spheres.
  virtual Sphere_Container *sphere_container() = 0;
  virtual bool balance_container() {return false;} // if needed
  virtual Rebalance_Stats *rebalance_stats() {return & sphere_container()->_rebalance_stats;}
  
  virtual void add_sphere( size_t si ) { sphere_container()->add_sphere(si); }
  
  Search_Structure(size_t num_dim) :  _tt(num_dim), _is_global(false) {}

  virtual ~Search_Structure()
  {}
  
  // statistics about the searches performed
  Search_Stats _stats;

protected:
  
  Trimming_Tool _tt;
  
  // if true, then this is a global search that iterates over *all* spheres, and the sphere_container contains all the spheres in the domain.
  // useful for high dimensions
  bool _is_global;

};

#endif
