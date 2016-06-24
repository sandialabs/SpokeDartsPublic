// Spoke_Darts.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// mid-level algorithms using spoke-darts

#ifndef SPOKE_DARTS
#define SPOKE_DARTS

#include "Geometry_Tool.hpp"
#include "Trimming_Tool.hpp"
#include "Spheres.hpp"
#include "Domain.hpp"
#include "Spoke_Length.hpp"
#include "Piercing.hpp"
#include "Random.hpp"

class Sphere_Container;
class Search_Structure;
class Wheel_Stats;

class Spoke_Darts : public Trimming_Tool
{
public:
  Spoke_Darts( Spheres &spheres, Domain &domain ) 
  : Trimming_Tool(domain.num_dim()), _spheres(spheres), _domain(domain), spoke_darts_workspace(domain.num_dim()), _rng( &Random::random_instance )
  {}

  void set_dimension( size_t num_dim )
  { Geometry_Tool::set_num_dim(num_dim); _spheres._pt.set_num_dim(num_dim); _domain.set_dimension(num_dim); }
  
  // many methods use a random number generator.
  // set it here if you don't want the static one.
  // This does not own the rng, the caller must new and delete it.
  void set_rng( Random *rng ) {_rng = rng;}
  Random *get_rng() { return _rng; }
  
  //===========================================================================
  // Interface for dart-based algorithms
  // use this if you want all the segments hit by a line, not just one of them
  // ==========================================================================
  // pick a segment
  // that is uncovered by trim_spheres, 
  // that is a subsegment of A,B, and [-A,-B] starting from point c and going in direction u
  // Note instead of returning a negative segment, this flips the sign of u, and A,B >= 0
  // returns false if there was no such segment
  bool pick_segment_piercing(const double *c, double * u, const double A, const double B,
                             Sphere_Container *trim_spheres,
                          double & seg_start, double & seg_end,
                          bool full_segment_only = false, double min_length = 0.);


  // search a circle of the sphere to find an uncovered point
  // will change u
  // return true if successful
  bool wheel(const double *c, const double rc, double * u, double * dart,
             Sphere_Container *trim_spheres, const double radius_factor = 1.);


  bool pick_segment_crossing(const double *c, const double rc, double *u, double *v,
                             Sphere_Container *trim_spheres, const double radius_factor,
      double &theta_start, double &theta_end, const double min_length = 0.);

  
  
  // find an uncovered dart by one throw, true if one was found
  bool generate_dart(double *dart, double *u, const double A, const double * c, Search_Structure *anchor_nbr,
                     const bool do_wheels, bool *dart_is_from_wheel_ptr, Wheel_Stats *wheel_stats,
                     bool *covered_sphere, const double radius_factor = 1.);

  // return some distance along the spoke length, between the interval of the spoke length
  double sample_from_spoke(Spoke_Length &sl, double r,
                           double sample_start_min, double sample_start_max, double mid_1, double mid_2,
                           double top, double mid_frac );
  
private:
  Spheres &_spheres;
  Domain &_domain;
  Random *_rng;
  //
  
  // workspace
  class Spoke_Darts_Workspace : public Point_Tool
  {
  public:
    // data
    Piercings piercings_p, piercings_m;
    std::vector<double> segments_p, segments_m;
    
    Crossings crossings;
    std::vector<double> arcs;
    
    // points
    double *v, *w, *c_s, *c_sp, *sp_s;
    
    ~Spoke_Darts_Workspace()
    {
      delete_sphere(v);
      delete_sphere(w);
      delete_point(c_s);
      delete_point(c_sp);
      delete_point(sp_s);
    }
    
    Spoke_Darts_Workspace(size_t num_dim) : Point_Tool(num_dim)
    {
      v = new_sphere();
      w = new_sphere();
      c_s = new_point();
      c_sp = new_point();
      sp_s = new_point();
    }
  };
  Spoke_Darts_Workspace spoke_darts_workspace;
};


  
#endif
