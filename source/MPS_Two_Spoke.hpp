// MPS_Two_Spoke.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef MPS_TWO_SPOKE_HPP
#define MPS_TWO_SPOKE_HPP

#include "MPS.hpp"
#include "Spoke_Length.hpp"

class MPS_Two_Spoke : public MPS
{
public:
  // =========
  // INTERFACE
  // call this function to create an MPS-like distribution using spokes
  // =========

  MPS_Two_Spoke() : MPS() { _nominal_beta = 3.; _max_beta = 4.; _shell_distance = 5.2; _pack_distance = 8.; }

  void create(size_t num_dim, bool is_periodic, double r, Search_Factory::Search_Type search_type = Search_Factory::UNSPECIFIED, bool all_searches_global = false, bool perm_ghosts = false);

//  static double max_radius(bool is_periodic) { return (is_periodic) ? 0.071 : 1.;}
  static double max_radius(bool is_periodic) { return (is_periodic) ? 1. : 10.;} // 1000 is a proxy for infinity

protected:

  // return some distance along the spoke length, between the interval of the spoke length
  double sample_from_spoke(Spoke_Length &sl, double r);


// for an uncovered spoke (first spoke), where to put the free spoke center (second spoke),
  // as a multiple of 2r
  const double _covered_mid_1 = 1.40;
  const double _covered_mid_2 = 1.70;
  const double _covered_top = 2.0;

  // for a free spoke (second spoke), where to put the sample point (covering disk center)
  // as a multiple of r
  const double _free_mid_1 = 1.30;
  const double _free_mid_2 = 1.70;
  const double _free_top = 2.0;
  
  // midpoint for a closed spoke
  const double _mid_frac = 0.4; // 0.4

  // max number of free spokes to throw, picking one of them for the new sample point
  static const size_t _num_free_spokes = 5;

  virtual std::string standard_prefix( size_t dim, bool periodic, double r, size_t loop_count )
  { return MPS::standard_prefix(dim, periodic, r, loop_count, "2two"); }

};

#endif