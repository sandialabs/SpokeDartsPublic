// Search_Range.hpp
// The search routines that use range trees

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0


#ifndef SEARCH_RANGE_HPP
#define SEARCH_RANGE_HPP

#include "Search_Structure.hpp"
#include "Range_Tree.hpp"
#include "Trimming_Tool.hpp"


// This is the geometric search methods performed over range trees

// Class structure is a littler different than the array-based one, Search_Array,
// in that this *is* a Range_Tree, a Sphere_Container; and Search_Array *has* a Sphere_Container.

class Search_Range : public Search_Structure, public Range_Tree 
{
public:

  Search_Range(Spheres &spheres) : Range_Tree(spheres), Search_Structure(spheres.num_dim()) {}

  // See Search_Structure interface class for descriptions of these functions
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() );

  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.);

  virtual void all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false );

  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
    const double radius_factor = 1. );

  virtual Sphere_Container *sphere_container() {return this;} // the Range_Tree
  
  virtual bool balance_container()
  { if (needs_rebalance()) { rebalance(); return true; } return false;}

private:
  // recursive functions called by public non-recursive function
  // much of the code is redundant and could be removed with templates, or with function pointers, or wrapping with an iterator
  void nearest_sphere_leaves( Tree_Node *node, size_t d );
  void nearest_sphere_recurse( One_Tree *tree );
  void all_near_spheres_recurse( One_Tree *tree );
  void add_all_leaves( Tree_Node *node );
  void add_leaf( size_t isphere );

  bool no_near_spheres_recurse( One_Tree *tree );
  bool no_near_leaves( Tree_Node *node );
  bool no_near_leaf( size_t isphere );

  
  void trim_line_anchored_recurse( One_Tree *tree );
  void trim_all_leaves( Tree_Node *node );
  void trim_leaf( size_t isphere );
  
  // any sphere whose dth center coordinate is less than x_lo or greater than x_hi
  // can not intersect the interval (defined by the trim statespace)
  void get_line_neighborhood( size_t d, double &x_lo, double &x_hi );

  // return the first node that is straddled by x_lo and x_hi
  // (prior nodes had both x_lo and x_hi on either the right or left)
  //  if update_closest, then we update x_lo and x_hi to only search for nodes that might be closer 
  //    than the current closest node we've found so far
  Tree_Node * find_split( One_Tree *tree, double &x_lo, double &x_hi, bool update_closest = false );  
  
  // variables that remain constant during the recursions for near points
  // put in the class for reduced memory allocation or stack size within the recursion
  // used by old_lo_hi, new_lo_hi, new_old_lo_hi, nearest_sphere, nearest_sphere_recurse
  class Near_Workspace
  {
  public:
    // data
    const double * _p;
    double _thresh, _thresh_squared;
    Sphere_Container *_sub_container;
    size_t _conflict_sphere;
    size_t _num_calls; // debug, search stats
    
    ~Near_Workspace() {}
    
  };
  Near_Workspace near_workspace;

  // update the lo and hi range of the near workspace based on old, new, or both spheres
  // uses near_workspace
  void old_lo_hi( size_t d, double& x_lo, double &x_hi );
  bool new_lo_hi( size_t &node, size_t d, double &x_lo, double &x_hi);
  bool new_old_lo_hi( size_t &node, size_t d, double &x_lo, double &x_hi);
  
  // Keep track of the state of a trimmed spoke within the recursive trimming function,
  // Encapsulated here to avoid explicitly passing the data to functions
  // used by trim_line_anchored, get_line_neighborhood, trim_leaf
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
      assert( fabs(*_A_1) >= _r ); // the segment must be outside the initial sphere
      assert( fabs(*_A_2) >= _r ); // ? * _radius_factor?
    }
  };
  Trim_Statespace trim_statespace;
};

#endif