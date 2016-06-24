// MPS_Favored.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "MPS_Favored.hpp"

#include "Ghost_Spheres.hpp"
#include "Ghost_Global_Container.hpp"
#include "Spoke_Length.hpp"
#include "Random.hpp"
#include "Search_Tree.hpp"
#include "Nested_Searches.hpp"
#include "Darts_IO.hpp"

bool MPS_Favored::spoke_is_long( const Spoke_Length &sl, double r )
{
  if (sl.A_2 - sl.A_1 < _short_skip_max * r)
  {
    if (sl.A_2 - sl.A_1 < _short_skip_min * r)
      return false;
    else
    {
      double st = _rng->generate_a_random_number();
      // zzyk. research question
      // Here st is cubed, regardless of dimension. It was chosen experimentally for d=2.
      // Perhaps the power should be dimension dependent, num_dim+1?
      return (sl.A_2 - sl.A_1) > (st * st * st * (_short_skip_max - _short_skip_min) + _short_skip_min) * r;
    }
  }
  return true;
}


void MPS_Favored::create(size_t num_dim_in, bool is_periodic, double r, Search_Factory::Search_Type search_type, bool all_searches_global, bool perm_ghosts)
{
  std::cout << std::endl << "Starting favored-spokes d:" << num_dim_in << " r:" << r << std::endl;

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

  // _do_central_beta = true;

  // spokes go from 1r to 3.8r
  const double anchor_dist = 1.;
  const double spoke_extent = 3.8;
  
  const bool do_wheels = false;

  const bool do_greedy_pass = true; // should the algorithm throw darts greedily?
  bool greedy_pass = false; // is the current set of dart throws greedy
  const size_t greedy_max_misses = 6;
  
  // num_successive_misses = 6; from orig
  _max_misses = 6; // 6 or 12?

  assert( greedy_max_misses <= _max_misses ); // change loop control if we want a larger value for greedy misses

  //===============
  // Setup
  //===============

  double max_search_distance = (spoke_extent + 1.1) * r; // .1 is a safety factor
  _influence_neighborhood = spoke_extent+1.5; // for reporting histogram
  truncate_influence_neighborhood( r );

  MPS::create(num_dim_in, is_periodic, r, max_search_distance, search_type, 0, perm_ghosts, all_searches_global);

  // Global_Container &all_spheres = *(ghosts ? new Ghost_Global_Container(*ghosts) : new Global_Container(*_spheres));

  // Search Structures

  Spoke_Length sl(anchor_dist, spoke_extent, r);

  bool dart_is_from_wheel = false;
  Wheel_Stats wheel_stats;
  bool covered_sphere(false);

  _nested_searches.clear();
  _nested_searches.reserve(3);

  _nested_searches._searches.push_back( _global_search );

  // subset of global space
  const double spoke_nbr_dist = (spoke_extent + 1.) * r;
  Search_Structure *spoke_nbr = Search_Factory::new_search( search_type, _spheres, all_searches_global, spoke_nbr_dist, 1., 0.);
  //  if (!all_searches_global) // always
  {
    _nested_searches._searches.push_back( spoke_nbr );
    _nested_searches._distances.push_back( spoke_nbr_dist );
  }

  // subset of spoke_nbr
  const double anchor_nbr_dist = (anchor_dist + 1.) * r;
  Search_Structure *anchor_nbr = Search_Factory::new_search( search_type, _spheres, all_searches_global, anchor_nbr_dist, 1., 0.);
  //  if (!all_searches_global) // always
  {
    _nested_searches._searches.push_back( anchor_nbr );
    _nested_searches._distances.push_back( anchor_nbr_dist );
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

  // first disk
  create_center_sphere( dart, r );
  Global_Container front_spheres(*_spheres);
  bool failed(false);

  // report on progess periodically
  std::cout << "Starting Favored Spokes! " << std::endl;
  size_t next_report(0);
  report(1, next_report);
  next_report = 1000 / num_dim();

 
  // for all front spheres
  for ( _front_disk = front_spheres.first(); !failed && _front_disk != front_spheres.bad_sphere_index(); _front_disk = front_spheres.next() )
  {

    _outer_loop_count = _loop_count;
    if (/* DISABLES CODE */ (0) && (_outer_loop_count == 61476))
    {
      std::cout << "outer, loop count " << _loop_count << " debug me!" << std::endl;
      std::cout << "_front_disk: " << _front_disk << " ";
      _spheres->_pt.print_sphere( (*_spheres)[_front_disk] );
      std::cout << std::endl;
      Point_Tool::_verification_level = 1; //0, 1, 2
    }

    
    // ====================================
    // advance the front
    // ====================================
    if (skip_front_disk())
      continue;

    advance_front();
    greedy_pass = do_greedy_pass;

    // quit if we've produced the requested number of layers around the central sphere
    if (reached_layer_limit())
      break;

    // throw spoke darts
    do 
    {

      ++_loop_count;

      // debug
      if (/* DISABLES CODE */ (0) && ( _loop_count == 61478 || (*_spheres)[0] == 0) )
      {
        std::cout << "inner, loop count " << _loop_count << " debug me!" << std::endl;
        Point_Tool::_verification_level = 1;
      }
      
      // ====================================
      // pick anchor point
      // ====================================

      sl.A = r;      
      bool valid_dart = _sd->generate_dart(dart, u, sl.A, _c, anchor_nbr,
                                    do_wheels, &dart_is_from_wheel, &wheel_stats, &covered_sphere );
            // if the whole sphere is covered, we know all future throws will fail, too
      if (covered_sphere)
        _quit_early = true;


      // =====================================================
      // Trim spoke through uncovered anchor
      // =====================================================
      if (valid_dart)
      {
        sl.reset(r);

        _sd->trim_by_domain( _c, u, sl.A_2, &_domain ); // does nothing if periodic
        assert(sl.A_2 >= sl.A_1);
        
        spoke_nbr->trim_line_anchored(_c, u, r, sl.A, sl.A_1, sl.A_2);

        // greedy?
        if (greedy_pass && sl.was_trimmed())
        {
          valid_dart = false;
        }
        
        // skip short spokes
        else 
          valid_dart = spoke_is_long( sl, r );
      }
      
      // ==========================================================
      // place sample point on non-empty trimmed dart
      // ==========================================================
      if (valid_dart)
      {
        // pick a sample point from the interval [A_1, A, A_2], but non-uniformly
        // store it in "dart"
        double dx = _sd->sample_from_spoke(sl, r, _sample_start_min, _sample_start_max, _mid_1, _mid_2, _top, _mid_frac);

        // double dx = _rng->random_uniform(sl.A_1, sl.A_2);
        _sd->axpy(dart, dx, u, _c);
      }
      
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
            io.plot_vertices_2d_allblack( io.loop_string( _loop_count, "B"), _spheres, &_domain, 1. );
          }
          if (draw_report_shaded_disks)
          {
            if (valid_dart)
              io.plot_vertices_2d_dart(  io.loop_string(_loop_count), _spheres, &_domain, 1., _front_disk, dart, &sl, u);
            else
              io.plot_vertices_2d_nodart( io.loop_string(_loop_count), _spheres, &_domain, 1., _front_disk );
          }
        }
      }
      
      // ===================================================================
      // add the dart to the pool of spheres, and subsequent local searches
      // ===================================================================
      if (valid_dart)
      {
        failed = failed || !create_new_sphere( dart );
      }

      // ====================================
      // miss
      // ====================================
      if (!valid_dart)
      {
        ++_misses;

        // greedy pass
        if ( greedy_pass && _misses >= greedy_max_misses )
        {
          greedy_pass = false;
          _misses = 0;
        }
      }

      // debug
//      if ( _ghosts->num_ghosts() > 1000 )
//      {
//        std::cout << _ghosts->num_ghosts() <<  " ghosts exist" << std::endl;
//      }
      report(_loop_count, next_report);

    } while ( _misses < _max_misses && !_quit_early ); // for dart throws

  }

  // debug plots
  if (/* DISABLES CODE */ (0))
  {
    Darts_IO io;
    if (draw_report_solid_disks)
    {
      io.plot_vertices_2d_allblack( io.loop_string(_loop_count, "B"), _spheres, &_domain, 1. );
    }
    if (draw_report_shaded_disks)
    {
      io.plot_vertices_2d_nodart( io.loop_string(_loop_count), _spheres, &_domain, 1., _front_disk );
    }
  }

  winnow_disks_for_output();
  
  output_results();

  //cleanup
  _sd->delete_sphere(dart);
  _sd->delete_sphere(u);
  
  Search_Factory::delete_search( anchor_nbr);
  Search_Factory::delete_search( spoke_nbr );
  _nested_searches.clear();

}