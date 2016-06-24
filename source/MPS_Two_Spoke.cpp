// MPS_Two_Spoke.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "MPS_Two_Spoke.hpp"

#include "Ghost_Spheres.hpp"
#include "Ghost_Global_Container.hpp"
#include "Spoke_Length.hpp"
#include "Random.hpp"
#include "Search_Tree.hpp"
#include "Nested_Searches.hpp"
#include "Darts_IO.hpp"

void MPS_Two_Spoke::create(size_t num_dim_in, bool is_periodic, double r, Search_Factory::Search_Type search_type, bool all_searches_global, bool perm_ghosts)
{
  std::cout << std::endl << "Starting two-spokes d:" << num_dim_in << " r:" << r << std::endl;

  /*
  find uncovered point
     spoke search, with 2r circles.
        pick an uncovered point as before
  find free point
     spoke search, centered at result of prior spoke, using r circles, and 2r extent in either direction
      perform several, to seek a far point that will cover this
     select sample from one of the uncovered segments, using swept area, ramp-function
 */

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

  const size_t num_free_spokes_used = 1; // how many do I actually throw

  const bool do_wheels = false;

  // default _max_misses of 12 is fine
  
  // _do_central_beta = true;

  //===============
  // Setup
  //===============

  bool dart_is_from_wheel = false;
  Wheel_Stats wheel_stats;
  bool covered_sphere(false);

  double anchor_nbr_dist = 4. * r;
  double trim_nbr_dist = 6. * r; // 6 = uses 2r coverage disks, to check a point up to 4r away
  double free_nbr_dist = 7. * r; // 7 = 4r uncovered point, plus 2r spoke length, plus 1r sphere radius

  // Initial and trimmed spoke lengths
  // r_free == r
  const double r_cover = 2. * r;
  Spoke_Length covered_sl(1., 2., r_cover); // 2r coverage disks, spoke goes from 2r to 4r
  Spoke_Length free_sl(-2., 2., r);       // 1r free disks
  Spoke_Length free_sl_onesided(0, 2., r); // 1r free disks
  Spoke_Length free_spokes[_num_free_spokes];
  for (size_t i = 0; i < _num_free_spokes; ++i)
    free_spokes[i].reset(-2., 2., r);
  // be sure trimming retains both sides of anchor for the free spokes

  const double max_search_distance = free_nbr_dist + .01 * r; // .01 is a safety factor

  _influence_neighborhood = (free_nbr_dist/r) + 1.;
  truncate_influence_neighborhood( r );
  
  MPS::create(num_dim_in, is_periodic, r, max_search_distance, search_type, 0, perm_ghosts, all_searches_global);

  // Global_Container &all_spheres = *(ghosts ? new Ghost_Global_Container(*ghosts) : new Global_Container(*_spheres));

  // Search Structures

  _nested_searches.clear();
  _nested_searches.reserve(4);

  _nested_searches._searches.push_back( _global_search );

  // subset of global space
  Search_Structure *free_nbr = Search_Factory::new_search( search_type, _spheres,  all_searches_global, free_nbr_dist, 1., 0.);
//  if (!all_searches_global) // zzyk ? do always
  {
    _nested_searches._searches.push_back( free_nbr );
    _nested_searches._distances.push_back( free_nbr_dist );
  // be sure to add the front_disk,  add_sphere( front_disk )
  }
  
  // subset free space
  Search_Structure *trim_nbr = Search_Factory::new_search( search_type, _spheres,  all_searches_global, trim_nbr_dist, 1., 0.);
  if (!all_searches_global)
  {
    _nested_searches._searches.push_back( trim_nbr );
    _nested_searches._distances.push_back( trim_nbr_dist );
  }

  // subset of spoke_nbr
  Search_Structure *anchor_nbr = Search_Factory::new_search( search_type, _spheres,  all_searches_global,  anchor_nbr_dist, 1., 0.);
  if (!all_searches_global)
  {
    _nested_searches._searches.push_back( anchor_nbr );
    _nested_searches._distances.push_back( anchor_nbr_dist );
  }

  double *dart = _sd->new_sphere();  // first dart
  double *dart2 = _sd->new_sphere();  // second dart
  double *u = _sd->new_sphere();  // first spoke direction
  double *u2 = _sd->new_sphere(); // second spoke direction
  double *u2m = _sd->new_sphere(); // minus second spoke direction, for some trimming 
  double *uncovered_point =  _sd->new_sphere();

  //===============
  // debug
  //===============
  _loop_count = 0;
  _outer_loop_count = 0;
  Point_Tool::_verification_level = 0;
  const bool do_debug = false;
  bool do_figures = (num_dim_in == 2);

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
    do_figures = false;
  }

  
  _setup_time.collect_stats();
  _main_time.start_clock();

  // first disk
  create_center_sphere( dart, r );
  _spheres->_pt.radius(dart2) = r; // fixed for all spheres currently
  Global_Container front_spheres(*_spheres);
  bool failed(false);

  // report on progess periodically
  std::cout << "Starting Two_Spoke MPS for step blue noise! " << std::endl;
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
    // beta-determining algorithm
    // just sample near the center of the sphere
    if (skip_front_disk())
      continue;
    
    advance_front();

    // quit if we've produced the requested number of layers around the central sphere
    if (reached_layer_limit())
      break;

    // the front_disk trims the free spoke, but is usually not in any of the search structures
    free_nbr->add_sphere(_front_disk);

    // throw spoke darts
    do 
    {

      ++_loop_count;
      
      const bool debug_me = false; // _loop_count == 4; // false; // (_loop_count >= 1562);
//      if (debug_me)
//        std::cout << "zzyk loop debug";
      
      // debug
      if ( (/* DISABLES CODE */ (0) && ( debug_me || (*_spheres)[0] == 0) ) )
      {
//        debug_me = true;
        std::cout << "inner, loop count " << _loop_count << " debug me!" << std::endl;
//        Point_Tool::_verification_level = 1; //0, 1, 2
      }
      
      // ====================================
      // pick anchor point, on 2r coverage disk
      // ====================================

      // dart on coverage sphere surface
      // variable radii
      // const double cr = _pt.radius(_c);
      // covered_sl.A = 2. * cr; 
      covered_sl.A = r_cover;  // constant radius disks
      
      // the final 2 parameter scales the r grabbed from the sphere coordinate, to be 2r coverage spheres
      bool valid_dart = _sd->generate_dart(dart, u, covered_sl.A, _c, anchor_nbr,
                                    do_wheels, &dart_is_from_wheel, &wheel_stats, &covered_sphere,
                                    2. );
      // if the whole sphere is covered, we know all future throws will fail, too
      if (covered_sphere)
        _quit_early = true;


      // =====================================================
      // Trim spoke through uncovered anchor
      // =====================================================
      if (valid_dart)
      {
        covered_sl.reset(r_cover);

        _sd->trim_by_domain( _c, u, covered_sl.A_2, &_domain ); // does nothing if periodic
        assert(covered_sl.A_2 >= covered_sl.A_1);
        
        trim_nbr->trim_line_anchored(_c, u, r_cover, covered_sl.A, covered_sl.A_1, covered_sl.A_2, 2.);
      }
      
      // ================================================================
      // place free spoke center on trimmed dart
      // ================================================================
      if (valid_dart)
      {
        // pick a sample point from the interval [A_1, A, A_2], but non-uniformly
        // store it in "dart"
        double dx = _sd->sample_from_spoke(covered_sl, r_cover, 1., 1.,
                                           _covered_mid_1, _covered_mid_2, _covered_top, _mid_frac);


        // double dx = _rng->random_uniform(sl.A_1, sl.A_2);
        _sd->axpy(uncovered_point, dx, u, _c);
      }
      
      // =====
      // debug plot disks
      // =====
      if (do_debug)
      {
        if ( (_loop_count > next_report) || debug_me )
        {
          if (do_figures)
          {
            Darts_IO io;
            if (draw_report_solid_disks)
            {
              io.plot_vertices_2d_allblack( io.loop_string( _loop_count, "B"), _spheres, &_domain, 1. );
            }
            if (draw_report_shaded_disks)
            {
              io.set_draw_text(true);
              if (valid_dart)
                io.plot_vertices_2d_dart(  io.loop_string( _loop_count ), _spheres, &_domain, 1., _front_disk, dart, &covered_sl, u, uncovered_point);
              else
                io.plot_vertices_2d_nodart( io.loop_string(_loop_count), _spheres, &_domain, 1., _front_disk );
            }
          }
        }
      }

      // ==============================
      // second spoke
      //   dart is an uncovered point. shoot spokes from it.
      //   find a free placement for a disk to cover it.
      // ==============================
      if (valid_dart)
      {
        // pick dart from uncovered point
        for ( size_t shots = 0; shots < num_free_spokes_used; ++shots )
        {
          _rng->sample_uniformly_from_unit_sphere_surface( u2 );

          // trim spoke. the spoke is of length 2r, but exclusion disks are just 1r
          // trim_spoke trims both the bottom and the top by the domain boundary
          free_sl.reset(r);
          free_sl.A = 0.;  // anchor is at the uncovered point, not the point at distance r

          // we pass in r, the radius of *all* the spheres, in order for the search kdd tree to work correctly
          // however, we can trim down to A = 0, less than r
          free_nbr->trim_line_anchored(uncovered_point, u2, r, free_sl.A, free_sl.A_1, free_sl.A_2);
          assert( free_sl.A_1 <= free_sl.A && free_sl.A <= free_sl.A_2);
          assert( free_sl.A_1 >= free_sl.A_1_init && free_sl.A_2 <= free_sl.A_2_init);
          
          // trim both ends of the spoke by the domain
          if (!is_periodic)
          {
            _sd->trim_by_domain( uncovered_point, u2, free_sl.A_2, &_domain ); // does nothing if periodic
            assert( free_sl.A_1 <= 0. );
            
            _sd->assign( u2m, u2 ); // *p = *q
            _sd->multiply( u2m, -1.); // u2m = -u2
            double A_1m = - free_sl.A_1;
            _sd->trim_by_domain( uncovered_point, u2m, A_1m, &_domain ); // does nothing if periodic
            free_sl.A_1 = -A_1m;
            assert( free_sl.A_1 <= 0. );
          }
          
          // todo, pick the longest shot spoke or something like that
        }
        
        // sample uniformly by area from the double-ended trimmed spoke
        // longer term, generate several spokes, and pick uniformly from the collection
        {
          const bool is_minus = _rng->random_spoke_side( &free_sl, &free_sl_onesided );
          double dx = _sd->sample_from_spoke(free_sl_onesided, r, 0., 0.,
                                             _free_mid_1, _free_mid_2, _free_top, _mid_frac);
          if (is_minus)
            dx = -dx;
          _sd->axpy( dart2, dx, u2, uncovered_point);
        }
        
        // =====
        // debug plot disks
        // =====
        if (do_debug)
        {
          if ( (_loop_count > next_report) || debug_me )
          {
            if (do_figures)
            {
              Darts_IO io;
              if (draw_report_solid_disks)
              {
                io.plot_vertices_2d_allblack( io.loop_string( _loop_count, "2B"), _spheres, &_domain, 1. );
              }
              if (draw_report_shaded_disks)
              {
                io.set_draw_text(true);
                if (valid_dart)
                  io.plot_vertices_2d_dart2(  io.loop_string( _loop_count, "2"), _spheres, &_domain, 1., _front_disk, dart, &covered_sl, u, uncovered_point, uncovered_point, &free_sl, u2, dart2);
                else
                  io.plot_vertices_2d_nodart( io.loop_string( _loop_count, "2"), _spheres, &_domain, 1., _front_disk );
              }
            }
          }
        }
      }

      // =============================================
      // add the two-radii dart to the pool of spheres
      // disk-free, and providing unique coverage
      // =============================================
      if (valid_dart)
      {
        failed = failed || !create_new_sphere( dart2 );
      }

      // ====================================
      // miss
      // ====================================
      if (!valid_dart)
      {
        ++_misses;
      }

      // debug
//      if ( _ghosts->num_ghosts() > 1000 )
//      {
//        std::cout << _ghosts->num_ghosts() <<  " ghosts exist" << std::endl;
//      }
      if( valid_dart )
        report(_loop_count, next_report);

    } while ( _misses < _max_misses && !_quit_early ); // for dart throws

  }

  // debug plots
  if (/* DISABLES CODE */ (0) && do_figures)
  {
    Darts_IO io;
    if (draw_report_solid_disks)
    {
      io.plot_vertices_2d_allblack( io.loop_string( _loop_count, "B"), _spheres, &_domain, 1. );
    }
    if (draw_report_shaded_disks)
    {
      io.plot_vertices_2d_nodart( io.loop_string(_loop_count), _spheres, &_domain, 1., _front_disk );
    }
  }

  winnow_disks_for_output();
  
  output_results();

  // cleanup
//  _sd->delete_sphere(st);
//  _sd->delete_sphere(end);
  _sd->delete_sphere(dart);
  _sd->delete_sphere(u);
  _sd->delete_sphere(uncovered_point);
  
  Search_Factory::delete_search( anchor_nbr);
  Search_Factory::delete_search( trim_nbr );
  Search_Factory::delete_search( free_nbr );
  _nested_searches.clear();

}