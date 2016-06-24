// Search_Grid.hpp
// The search routines that use grids

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef SEARCH_GRID_HPP
#define SEARCH_GRID_HPP

#include "Search_Structure.hpp"
#include "Grid_Container.hpp"
#include "Trimming_Tool.hpp"

// This is the geometric search methods performed using a background grid

class Search_Grid : public Search_Structure, public Grid_Container
{
public:

  Search_Grid(Spheres &spheres, double search_distance, double xmax, double xmin)
  : Grid_Container(spheres, search_distance, xmax, xmin),
  Search_Structure(spheres.num_dim()), _ws(spheres.num_dim()) {}

  Search_Grid(Spheres &spheres, double search_distance, Domain *domain)
  : Grid_Container(spheres, search_distance, domain->xmax()[0], domain->xmin()[0]),
  Search_Structure(spheres.num_dim()), _ws(spheres.num_dim()) {}


  // See Search_Structure interface class for descriptions of these functions
  virtual double nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold = std::numeric_limits<double>::max() );

  virtual bool no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor = 1., const double dist = 0.);

  virtual void all_near_spheres( Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius = false );

  virtual void trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2,
    const double radius_factor = 1. );

  virtual Sphere_Container *sphere_container() {return this;} // the Grid_Container
  
protected:
  // find the closest point of the grid cell ci to the query point p, return in q
  void closest_point( Cell_Index ci, const double *p, double *q );

private:

  // recursion over adjacent grid cells  
  void nearest_sphere_recurse( size_t d );
  bool no_near_spheres_recurse( size_t d );
  void all_near_spheres_recurse( size_t d );


  // variables that remain constant during the recursions for near points
  // put here for reduced memory allocation or stack size within the recursion
  class Near_Workspace
  {
  public:
    // data
    const double * _p;
    double _pair_dist_threshold;
    Sphere_Container *_sub_container;
    size_t _conflict_sphere;
    size_t _num_calls; // debug, search stats
    size_t _distances_checked;
    Cell_Index _current_cell;
    
    Near_Workspace(size_t num_dim) : _current_cell(new int[num_dim]) {}
    ~Near_Workspace() { delete [] _current_cell;}
  };
  Near_Workspace _ws;

  friend class Search_Tester;

};

#endif