// MPS_Bridson.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "MPS_Bridson.hpp"

#include "Ghost_Spheres.hpp"
#include "Ghost_Global_Container.hpp"
#include "Spoke_Length.hpp"
#include "Random.hpp"
#include "Search_Tree.hpp"
#include "Nested_Searches.hpp"
#include "Darts_IO.hpp"

void MPS_Bridson::create(size_t num_dim_in, bool is_periodic, double r, Search_Factory::Search_Type search_type, bool all_searches_global, bool perm_ghosts)
{
  std::cout << std::endl << "Starting Bridson AFP annular sampling d:" << num_dim_in << " r:" << r << std::endl;
  
  // run time
  _run_time.start_clock();
  _setup_time.start_clock();
  
  //===============
  // Parameters
  //===============

  // Passed in:
  // num_dim
  // is_periodic
  // r

  // sampling annulus goes from 1r to 2r
  double anchor_dist = 1.;
  double spoke_extent = 2.;

  _max_misses = 30;

  // _do_central_beta = true;

  //===============
  // Setup
  //===============

  double max_search_distance = (spoke_extent + 1.1) * r; // .1 is a safety factor
  _influence_neighborhood = spoke_extent + 1.;
  truncate_influence_neighborhood( r );
  
  // creates ghost or regular spheres, and global_search with underlying_search
  MPS::create(num_dim_in, is_periodic, r, max_search_distance, search_type, 0, perm_ghosts, all_searches_global);

  // Global_Container &all_spheres = *(ghosts ? new Ghost_Global_Container(*ghosts) : new Global_Container(*_spheres));

  // Search Structures

  Spoke_Length sl(anchor_dist, spoke_extent, r);

  _nested_searches.clear();
  _nested_searches.reserve(1);
  _nested_searches._searches.push_back( _global_search );

  const double spoke_nbr_dist = (spoke_extent + 1.) * r;
  Search_Structure *spoke_nbr = Search_Factory::new_search( search_type, _spheres, all_searches_global, spoke_nbr_dist, 1., 0.);
  
//  if (!all_searches_global)
  {
    // subset of global space
    _nested_searches._searches.push_back( spoke_nbr );
    _nested_searches._distances.push_back( spoke_nbr_dist );
  }

  double *dart = _sd->new_sphere();
  double *u = _sd->new_sphere();

  //===============
  // debug
  //===============
  _loop_count = 0;
  _outer_loop_count = 0;
  Point_Tool::_verification_level = 0;
  bool draw_report_solid_disks = false;
  bool draw_report_shaded_disks = false;
  
  if (/* DISABLES CODE */ (0))
  {
    draw_report_solid_disks = false;
    draw_report_shaded_disks = false;
    _draw_solid_disks = (num_dim_in == 2);
    _draw_shaded_disks = (num_dim_in == 2);
    _draw_solid_disks = false;
    _draw_shaded_disks = false;
    _do_beta = (num_dim_in < 7);
    _draw_histogram = true;
  }
  // silent settings, for do_central_beta
  if (1)
  {
    draw_report_solid_disks = false;
    draw_report_shaded_disks = false;
    _draw_solid_disks = false;
    _draw_shaded_disks = false;
    _do_beta = false;
    _draw_histogram = false;
  }
  
  _setup_time.collect_stats();
  _main_time.start_clock();

  // ============
  // create first disk
  // ============
  create_center_sphere( dart, r );
  Global_Container front_spheres(*_spheres);
  bool failed(false);

  // report on progess periodically
  size_t next_report(0);
  report(1, next_report);
  next_report = 1000 / num_dim();
  
  // ====================================
  // for all front spheres
  // ====================================
  // note, in the original Bridson, we throw thirty darts, then put the point back into the queue if there was a hit!
  // the points are visited in random order, not depth first search!
  // Does this matter?
  for ( _front_disk = front_spheres.first(); !failed && _front_disk != front_spheres.bad_sphere_index(); _front_disk = front_spheres.next() )
  {

    // ====
    // debug
    // ====
    _outer_loop_count = _loop_count;
    if (/* DISABLES CODE */ (0) && (_outer_loop_count == 93))
    {
      std::cout << "outer, loop count " << _loop_count << " debug me!" << std::endl;
      std::cout << "_front_disk: " << _front_disk << " ";
      _spheres->_pt.print_sphere( (*_spheres)[_front_disk] );
      std::cout << std::endl;
//      Point_Tool::_verification_level = 1; //0, 1, 2
    }

    
    // ====================================
    // advance the front, resets miss count
    // ====================================
    if (skip_front_disk())
      continue;

    advance_front();
    
    // quit if we've produced the requested number of layers around the central sphere
    if (reached_layer_limit())
      break;

    // throw spoke darts
    do {

      ++_loop_count;

      // ====
      // debug
      // ====
      if (/* DISABLES CODE */ (0) && ( _loop_count == 128 ) ) // || (*_spheres)[0] == 0) )
      {
        std::cout << "inner, loop count " << _loop_count << " debug me!" << std::endl;
        Point_Tool::_verification_level = 2;
      }
      
      // ====================================
      // pick random point from annulus
      // ====================================
      sl.A = _rng->random_by_volume(sl.A_1, sl.A_2);
      _rng->sample_uniformly_from_unit_sphere_surface( u );
      _spheres->_pt.axpy( dart, sl.A, u, _c );
      
      // =====================================================
      // Check if it is inside the domain
      // =====================================================
      bool valid_dart = _domain.in_domain(dart);
      
      // =====================================================
      // Check if it is covered by a sphere
      // =====================================================
      size_t near_sphere = _spheres->bad_sphere_index();
      valid_dart = valid_dart && spoke_nbr->no_near_spheres(dart, near_sphere);

      // =====
      // debug plot disks
      // =====
      if (/* DISABLES CODE */ (0))
      {
        if (_loop_count > next_report)
        {
          Darts_IO io;
          if (draw_report_solid_disks)
          {
            io.plot_vertices_2d_allblack( io.loop_string( _loop_count, "B" ), _spheres, &_domain, 1. );
          }
          if (draw_report_shaded_disks)
          {
            if (valid_dart)
              io.plot_vertices_2d_dart(  io.loop_string( _loop_count ), _spheres, &_domain, 1., _front_disk, dart, &sl, u);
            else
              io.plot_vertices_2d_nodart( io.loop_string( _loop_count ), _spheres, &_domain, 1., _front_disk );
          }
        }
      }
      
      // =====================================================
      // If uncovered, accept it as a new sphere
      // add the dart to the pool of spheres, and subsequent local searches
      // =====================================================
      if (valid_dart)
      {
        failed = failed || !create_new_sphere( dart );
      }
      
      // =====================================================
      // Else increment misses
      // =====================================================
      if (!valid_dart)
      {
        ++_misses;
      }

      // ====
      // debug
      // ====
      report(_loop_count, next_report);

    } while ( _misses < _max_misses && !failed);
  }

  winnow_disks_for_output();
  
  output_results();

  // cleanup
  _sd->delete_sphere(dart);
  _sd->delete_sphere(u);
  
  Search_Factory::delete_search( spoke_nbr );
  _nested_searches.clear();

}
