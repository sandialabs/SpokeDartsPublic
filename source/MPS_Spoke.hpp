// MPS_Spoke.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef MPS_SPOKE
#define MPS_SPOKE

#include "MPS.hpp"

class MPS_Spoke : public MPS
{
public:
  // =========
  // INTERFACE
  // call this function to create an MPS-like distribution using spokes
  // =========

  MPS_Spoke() : MPS() { _nominal_beta = 1.5; _max_beta = 2.; _shell_distance = 3.; _pack_distance = 4.; }

  void create(size_t num_dim, bool is_periodic, double r, Search_Factory::Search_Type search_type = Search_Factory::UNSPECIFIED, bool all_searches_global = false, bool perm_ghosts = false);

//  static double max_radius(bool is_periodic) { return (is_periodic) ? 0.161 : 1.;}
  static double max_radius(bool is_periodic) { return (is_periodic) ? 1. : 10.;} // 10 is a proxy for infinity

protected:

  virtual std::string standard_prefix( size_t dim, bool periodic, double r, size_t loop_count )
  { return MPS::standard_prefix(dim, periodic, r, loop_count, "0spoke"); }


};


#endif