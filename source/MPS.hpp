// MPS.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef MPS_HPP
#define MPS_HPP

#include <cmath>
#include <cstddef>
#include "Timing_Stats.hpp"
#include "Domain.hpp"
#include "Spoke_Darts.hpp"
#include "Nested_Searches.hpp"
#include "Random.hpp"

#include "Search_Factory.hpp"

class Spheres;
class Ghost_Spheres;
class Search_Structure;
class Global_Container;

// deriving from this is just a way to flag a top-level algorithm
// that creates sphere distributions
// also it provides some commonly used control variables


class MPS
{
public:
  
  MPS( size_t num_dimension = 2 ) :
  _front_disk(0), _misses(0), _max_misses(12),
  _quit_early(false), _max_num_spheres(0), _quit_spheres(0),
  _domain(num_dimension, false), _sd(0), _spheres(0), _ghosts(0), _global_search(0), _underlying_search(0),
  _next_rebalance(16),
  _draw_solid_disks(false), _draw_shaded_disks(false), _draw_histogram(true), _influence_neighborhood(4),
  _input_radius(0.), _nominal_beta(2), _max_beta(4),
  _shell_distance(3), _pack_distance(4),
  _neighbor_time(0.), _neighbor_size(0.), _neighbor2_time(0.), _neighbor2_size(0.), _all_searches_global(false),
  _do_plain(true), _do_beta(true), _do_performance(true), _do_central_beta(false), _skip_distance2(1000000000), _neighbor_radius(1000000000),
  _loop_count(), _outer_loop_count(), _layer_limit( (size_t) -1), _current_layer(0), _layer_boundary(0), _rng(&Random::random_instance)
  {}

  virtual ~MPS() 
  {
    clear();
  }
  
  virtual void clear()
  {
    if (_underlying_search != _global_search)
      delete _underlying_search;
    delete _global_search;
    delete _spheres; // virtual, will delete ghosts spheres
    delete _sd;
    _underlying_search = 0;
    _global_search = 0;
    _spheres = 0;
    _ghosts = 0;
    _sd = 0;
  }

  // derived classes that override create should call MPS::create early.
  // max_num_spheres will be estimated if passed in as zero
  // WARNING: the implementation currently ignores the perm_ghosts flag
  virtual void create(size_t num_dim, bool is_periodic, double r, double max_search_distance,
                      Search_Factory::Search_Type search_type, size_t max_num_spheres = 0, bool perm_ghosts = false, bool all_searches_global = false);

  virtual void set_dimension(size_t num_dim, bool is_periodic)
  { if (_rng) _rng->set_num_dim(num_dim); if (_sd) _sd->set_dimension(num_dim); _domain.set(num_dim, is_periodic);}
  
  // io related
  void report_time(std::ostream &out);
  void report_darts(std::ostream &out, size_t num_darts);
  void report_search_stats(std::ostream &out);
  
  
  Spheres *output_spheres() { return _spheres; }
  size_t output_size() { return _spheres->num_real_spheres(); } // ignore ghosts
  double output_time() { return _main_time.cpu_time; } // ignore pre and post processing time
  double output_radius() { return _input_radius; }
  double nominal_beta() { return _nominal_beta; }
  double max_beta() { return _max_beta; }

  // return the maximum radius over which the method will actually run
  // for periodic domains, there must be only one periodic copy of a point that is within the largest search distance
  static double max_radius(bool is_periodic) { return (is_periodic) ? 1. : 1000.;} // 1000 is a proxy for infinity

  
  double neighbor_time() const {return _neighbor_time;}
  double neighbor_size() const {return _neighbor_size;}
  // time taken by the non-broadest searches, and the number found on average
  // if there is more than one sub-search, return the time-weighted average of the number of neighbors
  double neighbor2_time() const {return _neighbor2_time;}
  double neighbor2_size() const {return _neighbor2_size;}

  // return the distance to the farthest Voronoi vertex for the given sphere
  double estimate_beta( size_t sphere_index );
  // prefix is used to make unique file names so that concurent runs don't overwrite each other.
  double calculate_beta_with_qhull( size_t sphere_index, std::string prefix );

  
protected:

  // set the maximum number of real spheres (only) we can handle
  virtual void set_size( size_t max_num_spheres )
  { _max_num_spheres = max_num_spheres; set_quit_spheres(); }

  // an upper bound guess at the number of spheres we need
  virtual size_t packing_size_estimate(double r, double frame_size = 0.);

  size_t num_dim() const {return _domain.num_dim();}

  // new_spheres will allocate memory based on packing size estimate, if the passed in max_num_spheres == 0
  Spheres * new_spheres( bool do_ghosts, double r, double max_search_distance, size_t max_num_spheres = 0 );   // called by create
  void new_global_search( const double max_search_distance, Search_Factory::Search_Type search_type );

  void create_center_sphere(double *p, double r);

  // call this when there is a new _front_disk
  void advance_front();

  Domain _domain;
  Spoke_Darts *_sd;

  // try to avoid exposing whether we have ghost or real spheres, beyond their initial creation
  Ghost_Spheres *_ghosts;
  Spheres *_spheres;
  Search_Structure *_global_search;
  Search_Structure *_underlying_search;
  //Global_Container *_all_spheres, just create as needed from ghosts
  
  Timing_Stats _run_time, _setup_time, _main_time, _post_time;
  size_t _max_num_spheres, _quit_spheres;  // these are real spheres only
  bool _quit_early;
  size_t _front_disk;
  size_t _max_misses;
  size_t _misses;
  size_t _next_rebalance;
  
  // debug
  size_t _loop_count, _outer_loop_count;

  // if >=0, then only throw darts from the first and subsequent disks
public:
  size_t _layer_limit;
protected:
  size_t _current_layer, _layer_boundary;
  bool reached_layer_limit()
  { return ( _current_layer > _layer_limit); }

  void set_quit_spheres() { _quit_spheres = _max_num_spheres - (long) ceil(pow(3, num_dim())) - 2; }

  const double *_c; // the front disk, shorthand for_spheres[_front_disk]
  Nested_Searches _nested_searches;

  // add c and periodic copies to the spheres, ghost spheres, the global search, and the nested search
  bool create_new_sphere( const double * dart );
  
public:
  // these global values are control parameters
  // for which output files to write
  bool _do_beta, _do_plain, _do_performance;
  bool _draw_solid_disks;
  bool _draw_shaded_disks;
  bool _draw_histogram;
  double _influence_neighborhood; // max distance to report in histogram

protected:
  // limit the neighborhood by some fraction of the box size
  bool truncate_influence_neighborhood( double r );
  bool dump_disks(std::string s);

  // report progress of the main loop of an MPS algorithm
  void report(size_t loop_count, size_t & next_report);
  // create a standard string for keeping output files apart
  
  virtual std::string standard_prefix( size_t dim, bool periodic, double r, size_t loop_count )
  { return standard_prefix(dim, periodic, r, loop_count, "X"); }
  std::string standard_prefix( size_t dim, bool periodic, double r, size_t loop_count, std::string method_string );
  
  // nominal beta is the maximum beta we expect to see in practice
  // max beta is the theoretical maximum that is probabilistically guaranteed
  double _input_radius, _nominal_beta, _max_beta;
  double _shell_distance, _pack_distance;
  
  bool _all_searches_global;

  friend class Search_Tester;
  
  // saved values for the main algorithm
  void count_neighbors();
  double _neighbor_time, _neighbor_size, _neighbor2_time, _neighbor2_size;
  
  // beware counting the time and number of neighbors for the histogram
  // time taken by broadest search, and the number found on average
  double current_neighbor_time() const;
  double current_neighbor_size() const;
  // time taken by the non-broadest searches, and the number found on average
  // if there is more than one sub-search, return the time-weighted average of the number of neighbors
  double current_neighbor2_time() const;
  double current_neighbor2_size() const;
  
  
  // for just creating the spheres near the center
public:
  bool _do_central_beta; // default is false
  double _skip_distance2; // (_input_radius * 4.)^2;
  double _neighbor_radius;
protected:
  // true if _do_central_beta and distance from front disk center to the first disk is greater than _skip_distance
  bool skip_front_disk();

  // get rid of the disks we don't want to output
  void winnow_disks_for_output();
  
  // write performance stats, ps files, and disk data files
  void output_results();
  
public:
  // compute and write the beta of the cetral point directly
  void report_central_beta(std::string prefix = std::string(), bool estimate_beta = false );

private:
  
public:
  // get/set rng. Global static by default
  void set_rng( Random *rng ) {_rng = rng;}
  Random *get_rng() {return _rng;}
  
protected:
  // could be static, or own rng to avoid a static one in a threaded environment
  Random *_rng;
  
  
  
};


#endif
