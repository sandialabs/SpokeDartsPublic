// Sphere_Container.hpp 

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// Abstract class
//
// base class defining the common interface for tree and array based containers, 
// that contain indices into the "Spheres"
//
// it is also an iterator over itself.
// To get more iterators, use Sphere_Container_Interator

#ifndef SPHERE_CONTAINER
#define SPHERE_CONTAINER

#include "Spheres.hpp"
#include "Sphere_Container_Iterator.hpp"
#include "Timing_Stats.hpp"

class Sphere_Container: public Sphere_Container_Iterator
{
public:
  Sphere_Container(Spheres &spheres) : Sphere_Container_Iterator(spheres) {}
  virtual ~Sphere_Container() {}

  // common interface functions
  
  // get rid of the contents
  virtual void clear() = 0;

  // add one sphere to the container
  virtual void add_sphere(size_t isphere) = 0;

  // number of elements in the container
  virtual size_t size() const = 0;

  //==============================================
  // debug methods that are generically implemented
  //==============================================
  // true if all distances bigger than the radius

  // all pairwise, quadratic
  bool verify_min_distances( std::ostream &out = std::cout ) const;

  // all distances to p, linear
  bool verify_min_distances( const size_t pi, std::ostream &out = std::cout ) const
  { return verify_min_distances( _spheres[pi], pi, out ); }
  bool verify_min_distances( const double *p, const size_t pi = Spheres::bad_sphere_index(), std::ostream &out = std::cout ) const;

  // just between p and q
  bool verify_min_distances(const size_t pi, const size_t qi,
                            std::ostream &out = std::cout ) const
  { return verify_min_distances( _spheres[pi], _spheres[qi], pi, qi, out ); }
  bool verify_min_distances( const double *p, const double *q,
                            const size_t pi = Spheres::bad_sphere_index(),
                            const size_t qi = Spheres::bad_sphere_index(),
                            std::ostream &out = std::cout ) const;

  // list all the spheres in the structure. Note that only one virtual function needs to be overridden
  virtual void print( std::string name = std::string(), std::ostream &out = std::cout ) const; // override me
  virtual void print_distances( const double *c, std::string name = std::string(), std::ostream &out = std::cout, size_t max_lines = 0 ) const;


  // some implementations, such as trees, may need to be rebalanced sometimes
  virtual bool needs_rebalance() const { return false; }
  virtual void rebalance() {}
  Rebalance_Stats _rebalance_stats;

protected:

  bool verify_min_distances_loc( const double *p, const double *q,
                                size_t pi, size_t qi, // indices of p and q in _spheres, or bad index
      double &r, double&r_tol, double &r_squared, double &r_squared_tol,
      std::ostream &out = std::cout ) const;
  void get_min_distances( double &r, double&r_tol, double &r_squared, double &r_squared_tol) const;

};

#endif