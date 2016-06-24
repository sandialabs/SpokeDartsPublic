// MPS_Favored.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef MPS_FAVORED_HPP
#define MPS_FAVORED_HPP

#include "MPS.hpp"
#include "Spoke_Length.hpp"

class MPS_Favored : public MPS
{
public:
  // =========
  // INTERFACE
  // call this function to create an MPS-like distribution using spokes
  // =========

  MPS_Favored() : MPS() { _nominal_beta = 2.5; _max_beta = 4; /*not sure*/ _shell_distance = 3.8; _pack_distance = 4.4; }

  void create(size_t num_dim, bool is_periodic, double r, Search_Factory::Search_Type search_type = Search_Factory::UNSPECIFIED, bool all_searches_global = false, bool perm_ghosts = false);

//  static double max_radius(bool is_periodic) { return (is_periodic) ? 0.102 : 1.;}
  static double max_radius(bool is_periodic) { return (is_periodic) ? 1. : 10.;} // 1000 is a proxy for infinity


protected:

  // true if sl is longer than short_skip_max, or some random value between short_skip_min and short_skip_max
  bool spoke_is_long( const Spoke_Length &sl, double r );

  // return some distance along the spoke length, between the interval of the spoke length
  double sample_from_spoke(Spoke_Length &sl, double r, double sample_start_min, double sample_start_max);

  // determining whether a spoke is long enough
  // control parameters that are unchanging
  const double _short_skip_min = 2.2; // 2.0; //  1.8; // 2.24
  const double _short_skip_max = 2.5; //2.4; //  2.0; // 2.54 is huge, since the whole spoke is 3.8 - 1 = 2.8 long

  // placing a sample on a spoke  
  // max placement for an open spoke
  // _bot, random between these two extremes
  const double _sample_start_min = 1.3; // 1.5
  const double _sample_start_max = 1.8; // 1.8

  const double _mid_1 = 2.30; // 2.30
  const double _mid_2 = 2.70;  // 2.66
  const double _top = 3.4; // 3.4


  // midpoint for a closed spoke
  const double _mid_frac = 0.4; // 0.4    


  virtual std::string standard_prefix( size_t dim, bool periodic, double r, size_t loop_count )
  { return MPS::standard_prefix(dim, periodic, r, loop_count, "1favored"); }


};


#endif