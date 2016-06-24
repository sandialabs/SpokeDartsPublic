// Search_Tree.hpp
// The search routines that use KD_Trees

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0


#ifndef SEARCH_TREE_HPP
#define SEARCH_TREE_HPP

#include "Search_Structure.hpp"
#include "KD_Tree.hpp"
#include "Trimming_Tool.hpp"

// This is the geometric search methods performed over trees

// Class structure is a littler different than the array-based one, Search_Array,
// in that this *is* a KD_Tree, a Sphere_Container, ; and Search_Array *has* a Sphere_Container.

class Search_Tree : public Search_Structure, public KD_Tree 
{
public:

  Search_Tree(Spheres &spheres) : KD_Tree(spheres), Search_Structure(spheres.num_dim()) {}

  // See Search_Structure interface class for descriptions of these functions
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() );

  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.);

  virtual void all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false );

  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
    const double radius_factor = 1. );

  virtual Sphere_Container *sphere_container() {return this;} // the KD_Tree
  
  virtual bool balance_container()
  { if (needs_rebalance()) { rebalance(); return true; } return false;}
  
  virtual void add_sphere(size_t isphere) { KD_Tree::add_sphere(isphere); }

private:
  // recursive functions called by public non-recursive function
  void nearest_sphere_recurse( size_t d_index, size_t node );

  void all_near_spheres_recurse( size_t d_index, size_t node );

  bool no_near_spheres_recurse( size_t d_index, size_t node );

  void trim_line_anchored_recurse(size_t d_index, size_t sphere_index);

  void trim_line_all_recurse(size_t d_index, size_t sphere_index);

  
  // state of a trimmed spoke within the trimming algorithm.
  // encapsulated here to avoid having to pass this data within recursive functions.
  // used by trim_line_anchored and trim_line_anchored_recurse
  class Trim_Statespace
  {
  public:
    // data
    const double* _c;
    const double* _u;
    double _r;
    double _A;
    double *_A_1;
    double *_A_2;
    double _radius_factor;
    size_t _verification_level;
    size_t _visited, _found;
    
    //
    void set( const double* c, const double* u, const double r, const double A, double &A_1, double &A_2, const double radius_factor)
    { _c = c; _u = u; _r = r; _A = A; _A_1 = &A_1; _A_2 = &A_2; _radius_factor = radius_factor;}
    
    void anchored_trim(Trimming_Tool &tt, const double *s)
    {
      ++_visited;
      if ( tt.anchored_trim( _c, _u, s, _A, *_A_1, *_A_2, _radius_factor ) )
        ++_found;
      assert( *_A_1 <= *_A_2 ); // but either could be negative
      // This is only a valid check if we are trimming a spoke that started from a sphere
      //    assert( fabs(*_A_1) >= _r ); // the segment must be outside the initial sphere
      //    assert( fabs(*_A_2) >= _r ); // ? * _radius_factor?
    }
    
  };
  Trim_Statespace trim_statespace;

};

#endif