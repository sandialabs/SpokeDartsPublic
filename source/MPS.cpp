// MPS.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "MPS.hpp"

#include <stdlib.h>

#include "Spheres.hpp"
#include "Ghost_Spheres.hpp"
#include "Search_Structure.hpp"
#include "Search_Tree.hpp"
#include "Search_Array.hpp"
#include "Search_Grid.hpp"
#include "Search_Ghost.hpp"
#include "Search_Ghost_Array.hpp"
#include "Global_Container.hpp"
#include "Darts_IO.hpp"
#include "Quality_Tool.hpp"

void MPS::create(size_t num_dim, bool is_periodic, double r, double max_search_distance,
                         Search_Factory::Search_Type search_type, size_t max_num_spheres, bool perm_ghosts, bool all_searches_global)
{
  clear();
  
  _input_radius = r;
  
  // relevant if _do_central_beta
  const double neighbor_factor = _shell_distance;
  const double n_distance = _input_radius * neighbor_factor;
  _neighbor_radius = n_distance;
  const double skip_distance = _input_radius * _pack_distance;
  _skip_distance2 = skip_distance * skip_distance;

  
  // this is OK for arrays
  if (is_periodic && max_search_distance >= 1.0)
  {
    if (search_type == Search_Factory::UNSPECIFIED)
    {
      all_searches_global = true;
      search_type = Search_Factory::ARRAY;
    }
    
    if (search_type != Search_Factory::ARRAY || !all_searches_global)
    {
      // error, a sphere might overlap a periodic copy of itself
      std::cerr << "INPUT ERROR:" << std::endl;
      std::cerr << "Input radius " << _input_radius << " causes a search distance greater than the periodic domain size." << std::endl;
      std::cerr << "This is only supported for Search_Factory::ARRAY" << std::endl;
      std::cerr << "Rerun using a smaller r, or an aperiodic domain, or search structure Search_Factory::ARRAY and all_searches_global" << std::endl;
      assert(0);
    }
  }
  _all_searches_global = all_searches_global;
  set_dimension(num_dim, is_periodic);
  new_spheres( is_periodic, r, max_search_distance, max_num_spheres);
  new_global_search( max_search_distance, search_type );
}


size_t MPS::packing_size_estimate(double r, double frame_size)
{
  const double box_size = _sd->box_volume( _domain.xmax(), _domain.xmin(), frame_size );
  const double half_sphere_size = _sd->absolute_volume( r/2 ); // non-overlapping
  const double est_double  = box_size / half_sphere_size ;
  
  // account for at least one spheres per corner
  const double one_per_corner = pow( 2., (double) num_dim() );
  double upper_bound = est_double + one_per_corner;
  
  const size_t max_size_t = std::numeric_limits<size_t>::max() - 2;
  const size_t est = upper_bound > max_size_t ? max_size_t : (size_t) upper_bound;
  return est;
}


void MPS::create_center_sphere(double *p, double r)
{
  // coordinates = center of (periodic) domain
  for (size_t d = 0; d < num_dim(); d++)
  {
    // center
    p[d] = (_domain.xmin()[d] + _domain.xmax()[d] / 2.);
    // corner
//    p[d] = _domain.xmin()[d];
  }
  _sd->radius(p) = r;
  _global_search->add_sphere( _spheres->add_sphere(p) );
}

Spheres * MPS::new_spheres( bool do_ghosts, double r, double max_search_distance, size_t max_num_spheres )
{
  // create real or ghost spheres
  _ghosts = 0;
  bool perm_ghosts = false;
  if (do_ghosts)
  {
    // currently, it appears that temporary ghosts are always faster and take less memory
    // this might change if we successfully balance the search tree, or are using exhaustive search in very high dimensions
    perm_ghosts = num_dim() < 5;
    perm_ghosts = false;
    delete _ghosts;
    _ghosts = new Ghost_Spheres( _domain );
    _ghosts->set_frame( max_search_distance );
    _ghosts->set_permanent_ghosts( perm_ghosts );
    _spheres = _ghosts;
  }
  else
  {
    delete _spheres;
    _spheres = new Spheres(_domain);
  }
  _domain._bounds_tol = r * .1e-6;
  
  // spoke_darts over these spheres
  delete _sd;
  _sd = new Spoke_Darts(*_spheres,_domain);
  _sd->set_rng( _rng );
  
  // memory for actual disks
  // space for how many disks, 1,000,000 ?
  size_t array_size = max_num_spheres;
  size_t real_guess = max_num_spheres / 2;
  if (max_num_spheres == 0)
  {
    real_guess = packing_size_estimate(r);
    // plus one for each corner, plus center sphere
    size_t bonus = ( size_t(1)  << _domain.num_dim() ) + 1;
    array_size = real_guess + bonus;
    if (do_ghosts)
    {
      array_size = ( perm_ghosts ) ? packing_size_estimate(r, max_search_distance) : real_guess * 2;
    }
  }
  // 10 billion memory limit
  const size_t mem_limit = 10000000000;
  if (array_size > mem_limit)
  {
    array_size = mem_limit;
    std::cout << "Warning, memory limit of " << mem_limit << " spheres may be inadequate." << std::endl;
  }
  if (_do_central_beta && array_size > 200000)
    array_size = 200000;
  _spheres->new_sphere_array(array_size);
  set_size(array_size);
  
  // give the user some hint of how long this is going to take
  // time = O(2^d n d), where
  // n is the number of output spheres
  // 2^d is roughly proportional to N the number of neighbors. N is also bounded by n.
  // d is the time for every geometric operation such as taking distances, projections, storing a point, etc
  std::cout << "Creating MPS in dimension " << num_dim() << " with radius " << r << std::endl;
  std::cout << "Time and memory guess" << std::endl;
  std::cout << "Estimating at most " << real_guess << " spheres";
  if (do_ghosts)
  {
    std::cout << " plus " << array_size - real_guess << (perm_ghosts ? " permanent" : " temporary" ) << " ghosts.";
  }
  const double neighbor_guess =  pow(1.8, num_dim()-1);
  const double time_factor = 1.4e-5;
  // for some reason in practice it takes one less factor of d
  const double time_low = time_factor * num_dim() * neighbor_guess * array_size;
  const double time_high = time_factor * num_dim() * array_size * array_size;
  // std::cout << std::endl << "Relative time between " <<  time_low << " and " << time_high << ", closer to the second for higher dimensions." << std::endl;
  std::cout << std::endl << "Relative time estimate is " <<  time_low << std::endl;
  
  return _spheres;
}

void MPS::new_global_search( const double max_search_distance, Search_Factory::Search_Type search_type )
{
  assert( _spheres );

  // pick a method well suited to the dimension and search distance
  // if a grid, how many cells would we have?
  double num_cells(0);
  const double max_cells = 1.0e7;
  const double domain_side = _domain.xmax()[0] - _domain.xmin()[0]; // assumes a cube, not a rectangle
  const double domain_diagonal = sqrt(num_dim()) * domain_side;
  // heuristic choice, make cell side half way between the domain diagonal and domain side, times the search diameter
  double cell_side = 2. * max_search_distance; //  * ( domain_diagonal + domain_side ) / ( 2. * domain_side );
  num_cells = ceil( pow( domain_side / cell_side, num_dim()) );
  
  if (search_type == Search_Factory::UNSPECIFIED)
  {
    if ( max_search_distance > 0.3 * domain_diagonal )
    {
      Search_Array *search_array = new Search_Array(*_spheres, true);
      _underlying_search = search_array;
      search_type = Search_Factory::ARRAY;
    }
    else if ( num_cells < max_cells )
    {
      Search_Grid *grid = new Search_Grid( *_spheres, cell_side, _domain.xmax()[0], _domain.xmin()[0] );
      _underlying_search = grid;
      search_type = Search_Factory::GRID;
    }
    else
    {
      Search_Tree *tree = new Search_Tree( *_spheres );
      _underlying_search = tree;
      search_type = Search_Factory::TREE;
    }
  }
  else
    _underlying_search = Search_Factory::new_search(search_type, _spheres, true, cell_side, _domain.xmax()[0], _domain.xmin()[0] );
  
  _global_search = _underlying_search;
  
  // wrap the search with ghosts?
  if (_ghosts)
  {
    // max-search = half the domain width usually
    Search_Ghost_Structure * sg = 0;
    if (search_type == Search_Factory::ARRAY && _ghosts->is_temporary_ghosts() )
    {
      // fundamentally temporary ghosts
      Search_Array *search_array = dynamic_cast<Search_Array*>(_underlying_search);
      Search_Ghost_Array *sg_array = dynamic_cast<Search_Ghost_Array*>(_underlying_search);
      if ( search_array )
        sg = new Search_Ghost_Array( *_ghosts, search_array );
      else if ( sg_array )
        sg = sg_array;
      else
        assert(0); // unknown
    }
    else
    {
      sg = new Search_Ghost( *_ghosts, _underlying_search );
      // if ghosts are temporary, then auto discard and regather ghosts on each global search
      sg->set_search_reghost( true );

      const double max_OK_search = domain_side / 2.;
      if (  max_search_distance >= max_OK_search )
      {
        // input error, report a message to the user
        const double f(max_search_distance / max_OK_search);
        const double max_r( _input_radius / f);
        std::cerr << "INPUT ERROR:" << std::endl;
        std::cerr << "Radius is too big for ghosting to work, by a factor of " << f << std::endl;
        std::cerr << "Please reduce radius from " << _input_radius << " to at most " << max_r << std::endl;
        std::cerr << "Or make the domain aperiodic" << std::endl;
        std::cerr << "Or use search_ghost_arrays for your proximity searches in trimming" << std::endl;
        assert( max_search_distance < max_OK_search ); // always fails
      }

    }
    _global_search = sg;
  }
}

void MPS::advance_front()
{
  _c = (*_spheres)[_front_disk]; // shorthand
  
  if (_front_disk == 0)
  {
    _current_layer = 0;
    _layer_boundary = 1;
  }
  else if (_front_disk == _layer_boundary)
  {
    ++_current_layer;
    _layer_boundary = _spheres->num_real_spheres();
  }

  // rebalance, if needed, but not too often.
  if ( (_front_disk > _next_rebalance) && _global_search->balance_container() )
  {
    _next_rebalance = (size_t) (_front_disk * 2); // try setting the threshold based on _max_depth increasing by a factor of 2 ? 
  }
  
  if (!_all_searches_global)
    _nested_searches.gather_neighbors( _c );
  _misses = 0;
  _quit_early = false;


  // debug list the neighbors
  if (Point_Tool::_verification_level>=2)
    for ( size_t i = (_ghosts ? 1 : 0); i < _nested_searches._searches.size(); ++i )
    {
      std::cout << std::endl << "search level " << i << std::endl;
      _nested_searches._searches[i]->sphere_container()->print();
      std::cout << std::endl;
    }

  
  //== debug
  if (Point_Tool::_verification_level >= 1)
  {
    for ( size_t i = (_ghosts ? 1 : 0); i < _nested_searches._searches.size(); ++i )
    {
      if (!_nested_searches._searches[i]->sphere_container()->verify_min_distances( _front_disk ))
      {
          std::cout << "advance_front: distances to _front_disk " << _front_disk << " violated by search structure " << i << std::endl;
          assert(0);
      }
    }
  }
  if (Point_Tool::_verification_level >= 2)
  {
    if (_ghosts && !_nested_searches._searches[0]->sphere_container()->verify_min_distances( _front_disk ))
    {
      std::cout << "advance_front: distances to _front_disk " << _front_disk << " violated by global search structure " << 0 << std::endl;
      assert(0);
    }
    for ( size_t i = (_ghosts ? 1 : 0); i < _nested_searches._searches.size(); ++i )
    {
      if (!_nested_searches._searches[i]->sphere_container()->verify_min_distances())
      {
        std::cout << "advance_front: distances to two neighbor disks of " << _front_disk << " violated by search structure " << i << std::endl;
        assert(0);
      }
    }
  }
  //== end debug
}


bool MPS::create_new_sphere( const double * dart )
{
  size_t local_spherei, real_spherei;
  local_spherei = _spheres->add_sphere(dart, real_spherei);
  
  // the first search is a global one, and needs to know about *all* the new spheres
  // the other searches are local, and only need to know about the new local sphere
  const size_t first_local_search = _ghosts ? 1 : 0;
  if ( first_local_search )
  {
    if ( _ghosts->is_temporary_ghosts() )
    {
      _global_search->add_sphere( real_spherei );
      // any ghost isn't important, since the next global search will create the ghosts we need

      //== debug
      if (Point_Tool::_verification_level >= 1)
      {
        if (!_global_search->sphere_container()->verify_min_distances( real_spherei ))
        {
          std::cout << "distances to new real point " << real_spherei << " violated by global search structure " << std::endl;
          assert(0);
        }
      }
      //== end debug
    
    }
    else
      for ( size_t j = 0; j < _ghosts->_num_added_spheres; ++j )
      {
        const size_t sj = _ghosts->_added_spheres[j];
        _global_search->add_sphere( sj );
        
        //== debug
        if (Point_Tool::_verification_level >= 1)
        {
          if (!_global_search->sphere_container()->verify_min_distances( sj ))
          {
            std::cout << "distances to new point " << sj << " violated by global search structure" << std::endl;
            assert(0);
          }
        }
        //== end debug

      }
    
    //== debug
    if (Point_Tool::_verification_level >= 2)
    {
      if (!_global_search->sphere_container()->verify_min_distances())
      {
        std::cout << "distances between points violated by global search structure" << std::endl;
        assert(0);
      }
    }
    //== end debug


  }
  for ( size_t i = first_local_search; i < _nested_searches._searches.size(); ++i )
  {
    _nested_searches._searches[i]->add_sphere(local_spherei);

    
    //== debug
    if (Point_Tool::_verification_level >= 1)
    {
      if (!_nested_searches._searches[i]->sphere_container()->verify_min_distances( local_spherei ))
      {
        std::cout << "distances to new point " << local_spherei << " violated by search structure " << i << std::endl;
        assert(0);
      }
    }
    if (Point_Tool::_verification_level >= 2)
    {
      if (!_nested_searches._searches[i]->sphere_container()->verify_min_distances())
      {
        std::cout << "distances between points violated by search structure " << i << "" << std::endl;
        assert(0);
      }
    }
    //== end debug

  }
  
  // reset misses, and check for quitting
  _misses = 0;
  
  if ( _spheres->size() > _quit_spheres )
  {
    std::cout << "ERROR: failed, maximum number of spheres reached = " << _spheres->size() << ". Rerun with a larger max_num_spheres." << std::endl;
    _quit_early = true;
    return false; // i failed
  }
  return true;
}


bool MPS::dump_disks(std::string s)
{
  bool error(false);
  
  std::stringstream ss;
  
  // final arrangment
  Darts_IO io;
  if ( num_dim() > 1)
  {
    if (_draw_shaded_disks)
    {
      ss.str(std::string()); ss.clear(); ss << s << "-rings";
      io.plot_vertices_2d_rings(ss.str(), _spheres, &_domain, 1. );
    }
    
    if (_draw_solid_disks)
    {
      // should be all black
      ss.str(std::string()); ss.clear(); ss << s << "-solid";
      io.plot_vertices_2d_allblack(ss.str(), _spheres, &_domain, 1. );
      // should have no overlaps
      ss.str(std::string()); ss.clear(); ss << s << "-half";
      io.plot_vertices_2d_allblack(ss.str(), _spheres, &_domain, 0.5 );
      // should look interesting
      ss.str(std::string()); ss.clear(); ss << s << "-dots";
      io.plot_vertices_2d_allblack( ss.str(), _spheres, &_domain, 0.1 );
    }
  }
  
  if (_draw_histogram)
  {
    Quality_Tool quality_tool;
    Histogram histogram;
    error = quality_tool.build_histogram(histogram, _spheres, &_domain, _global_search, _influence_neighborhood) || error;
    io.RDF_histogram(std::cout, histogram );
    io.RDF_histogram(s, histogram);
    io.RDF_histogram_data(s, histogram);
  }

  
  return error;
}

double MPS::current_neighbor_time() const
{
  return _nested_searches._searches[0]->_stats.cpu_time;
}
double MPS::current_neighbor_size() const
{
  return _nested_searches._searches[0]->_stats.average_neighbors();
}

double MPS::current_neighbor2_time() const
{
  double t = 0.;
  for (size_t i = 1; i < _nested_searches._searches.size(); ++i)
    t += _nested_searches._searches[i]->_stats.cpu_time;
  return t;
}
double MPS::current_neighbor2_size() const
{
  const double total_t = neighbor2_time();
  double t(0.), N(0.), Nave(0.);
  for (size_t i = 1; i < _nested_searches._searches.size(); ++i)
  {
    t = _nested_searches._searches[i]->_stats.cpu_time;
    N = _nested_searches._searches[i]->_stats.average_neighbors();
    Nave +=  N * ( t / total_t);
  }
  return Nave;
}


void MPS::count_neighbors()
{
  _neighbor_time = current_neighbor_time();
  _neighbor_size = current_neighbor_size();
  _neighbor2_time = current_neighbor2_time();
  _neighbor2_size = current_neighbor2_size();
}


void MPS::report_time(std::ostream &out)
{
    out << "MPS run-time data:" << std::endl;
    out.precision(4);
    out<<  "    Total Time = " << std::fixed << _run_time.cpu_time << " seconds:" <<
      " setup = " << _setup_time.cpu_time << " & main = " << _main_time.cpu_time << " & post = " << _post_time.cpu_time << " )"
      << std::endl;
    //  handy format for spreadsheets
}

void MPS::report_darts(std::ostream &out, size_t num_darts)
{
  out<<  "    Number of created disks = " << std::fixed << _spheres->size();
  if (_domain.is_periodic())
    out << " = real:" << std::fixed << _spheres->num_real_spheres() << " ghosts:" << _ghosts->num_ghosts();
  out << std::endl;
  out<<  "    Number of thrown darts = " << std::fixed << num_darts << std::endl;

}

void MPS::report_search_stats(std::ostream &out)
{
  _nested_searches.report_stats(out);
}

void MPS::report(size_t loop_count, size_t & next_report)
{
  if (loop_count > next_report)
  {
    _main_time.collect_stats();
    _main_time.start_clock();
    std::cout << "elapsed " << _main_time.cpu_time << "sec, front sphere " << _front_disk << " of " << _spheres->size();
    if ( _ghosts && _ghosts->size() )
      std::cout << " (+ " << _ghosts->num_ghosts() << " ghosts )";
    std::cout << ", loop " << loop_count;
    std::cout << std::endl;
    next_report *= 1.5;
  }
}

//
//wheels
//    if (do_wheels)
//      out << " (spin time = " << wheel_stats.cpu_time << ", success time = " << wheel_stats.cpu_time_success << " )";
//    out << std::endl;
//    out.precision(0);
//    if (do_wheels)
//    {
//      out<<  "    Number of spun wheels = " << std::fixed << wheel_stats.spins << ", successes = " << std::fixed << wheel_stats.spins_success << std::endl;
//    }
//
//piercings
//    if (do_pierce)
//    {
//      out.precision(6);
//      const double f_success = ((double) pierce_stats.success) / ((double) pierce_stats.attempts);
//      const double f_u_success = ((double) pierce_stats.unique_success) / ((double) pierce_stats.success);
//      out << " Piercing rays, attempts: " << pierce_stats.attempts << ", fraction success: " << f_success << ", unique-success:" << f_u_success << std::endl;
//    }


bool MPS::truncate_influence_neighborhood( double r )
{
  // influence * rmust be less than 1/2 the box in order to function
  const double box_frac = 0.3;
  const double max_influence = box_frac / r;
  if (_influence_neighborhood > max_influence)
  {
    std::cout << "Warning, truncating " << _influence_neighborhood << "*r influence distance to " << max_influence << "*r. " << std::endl;
    _influence_neighborhood = max_influence;
    return true;
  }
  return false;
}


std::string MPS::standard_prefix( size_t dim, bool periodic, double r, size_t loop_count, std::string method_string )
{
  std::stringstream ss;
  int r_inverse = int( floor(1./r) );
  std::string periodic_string;
  if (periodic)
    periodic_string = "P";
  else
    periodic_string = "n";
    
  ss << dim << "_" << periodic_string << "_" << r_inverse << "_" << method_string << "_" << loop_count;
  return ss.str();

}

bool MPS::skip_front_disk()
{
  if (!_do_central_beta)
    return false;
  
  if (_front_disk == 0)
    return false;
  
  const double dist2 = _spheres->_pt.distance_squared( (*_spheres)[_front_disk], (*_spheres)[0] );

  if (dist2 > _skip_distance2)
    return true;
 
  return false;
}


void MPS::winnow_disks_for_output()
{
  if (_ghosts)
    _ghosts->dematerialize_temporary_ghosts();

  if (_do_central_beta)
  {
    // remove any sphere that is too far away from the central sphere
    if (!_spheres->size())
      return;
    
    // beware: neighor_radius could be different than skip_radius
    const double neighbor_radius_squared = _neighbor_radius * _neighbor_radius;
    Global_Container *gc = Ghost_Global_Container::new_global_container( *_spheres );
    const size_t bad_sphere_index = gc->bad_sphere_index();
    size_t i = 1;
    while ( i != bad_sphere_index )
    {
      if ( _spheres->_pt.distance_squared( (*_spheres)[i], (*_spheres)[0] ) > neighbor_radius_squared )
      {
        _spheres->remove_sphere(i);
        // set gc to have the state that the next returned value is i
        gc->set(i);
      }
      else
      {
        i = gc->next();
      }
    }
    delete gc;
  }
}


void MPS::output_results()
{
  
  // final statement of the same format as the report
  // also collects the final stats on the main loop
  std::cout << std::endl << "Finished! " << std::endl << "spheres = " << _spheres->size() << std::endl;
  size_t next_report = _loop_count-1;
  report(_loop_count, next_report);
  count_neighbors();
  
  
  _post_time.start_clock();
  
  Darts_IO io;
  std::string prefix = standard_prefix( num_dim(), _domain.is_periodic(), output_radius(), _loop_count );
  
  std::ostream *perf_out(0);
  if (_do_performance)
  {
    // write algorithm performance stats
    // write info about searches before they are used for the rdf plot
    std::stringstream ss;
    ss << prefix << "_performance_final.txt";
    perf_out = & (io.new_ostream(false, ss) );
    report_search_stats( std::cout );
    report_search_stats( *perf_out );
  }
  
  // write output files
  {
    if (_do_plain)
      io.save_spheres_plain( _spheres, prefix );
    //io.save_spheres_for_psa( _spheres, ss2.str() );
    //io.save_spheres_for_spectrum( _spheres, ss2.str() );
    
    // Output for Beta
    // set frame to beta*r
    // Because beta*r is the probabilistically guaranteed r_cover, and there could be a Vor vertex
    //   near the domain boundary defined by a sphere almost beta*r outside the domain.
    if (_do_beta)
    {
      const double write_frame = _ghosts ? nominal_beta() * output_radius() : 0.;
      io.save_spheres_for_beta( _spheres, write_frame, prefix );
    }
  }
  
  // dump output, ps files
  // do last since the histogram can take a long time
  std::stringstream ss2;
  ss2 << prefix << "_final";
  dump_disks(ss2.str());
  
  _post_time.collect_stats();
  _run_time.collect_stats();
  
  // rest of performance data
  if (_do_performance)
  {
    std::ostream &out = std::cout;
    report_time(std::cout);
    report_darts(out, _loop_count);
  }
  if (_do_performance)
  {
    report_time(*perf_out);
    report_darts(*perf_out, _loop_count);
    // time and output size
    *perf_out << std::endl << _main_time.cpu_time << "  " << _spheres->num_real_spheres() << std::endl;
    io.delete_ostream( false, *perf_out );
  }
}

void MPS::report_central_beta(std::string prefix, bool do_estimate)
{

  const double beta = ( do_estimate ? estimate_beta( 0 ) : calculate_beta_with_qhull( 0, prefix ) );

  // write output
  std::stringstream os;
  os << prefix << "_beta.txt";
  std::ofstream beta_file( os.str(), std::ios::app ); // append
  beta_file << beta << std::endl;
}

// todo, open and close beta_file once


double MPS::calculate_beta_with_qhull( size_t sphere_index, std::string prefix )
{

  // link to qhull
  // compute vor vertices (for cental point?)
  // find the vertices for the central point (only)
  // report the max distance encountered
  // append to file
  
  // for now, use a file interface, then see if it is fast enough
  
  assert(sphere_index == 0);
  
  // write beta_points.txt
  Darts_IO io;
  std::stringstream bs;
  bs << prefix << "_central";
  io.save_spheres_for_beta( _spheres, 0., bs.str() );
  
  std::stringstream vf;
  vf << prefix << "_voronoi.dat";
  
  // generate Voronoi vertices of all points
  // first source .bashrc to get the path to qvoronoi
  std::stringstream ss;
  ss << "source ~/.bashrc ; cat " << bs.str() << "_beta_points.txt | qvoronoi o TO \"" << vf.str() << "\"";
  /* int ok = */ system( ss.str().c_str() );
  //  /* int ok = */ system( "source ~/.bashrc ; cat central_beta_points.txt | qvoronoi o TO \"voronoi.dat\"" );
  
  // read in the Voronoi vertices file
  std::ifstream vor_file( vf.str(), std::ios::in );
  size_t vor_dim, num_vor_cells, num_vor_vertices, one;
  vor_file >> vor_dim >> num_vor_vertices >> num_vor_cells >> one;
  
  // number of voronoi cells should be equal to the number of disks
  assert( num_vor_cells == _spheres->size() );
  
  
  // read vor vertex coordinates
  // first one is the point at infinity
  Sphere_Array vor_vertices( num_dim(), num_vor_vertices );   // space for the coordinates
  for ( size_t i = 0; i < num_vor_vertices; ++i )
  {
    for ( size_t d = 0; d < num_dim(); ++d )
      vor_file >> vor_vertices[i][d];
  }
  
  
  // read first vor region, only!
  std::vector<size_t> cell0;
  size_t vi;
  std::string line;
  getline(vor_file, line);  // read spaces to endln, ready to read a real line
  
  // while(getline(vor_file, line)) // all lines
  getline(vor_file, line);  // first line only
  {
    std::istringstream vss(line);
    size_t cell_size;
    vss >> cell_size;
    while(vss >> vi)
    {
      cell0.push_back(vi);
    }
    assert( cell0.size() == cell_size );
  }
  
  // process the vor vertices to find beta
  const double ignore_d2 = _neighbor_radius * _neighbor_radius / 4.;
  double max_d2 = 0.;
  const double *c = (*_spheres)[0]; // center sphere
  for ( size_t v = 0; v < cell0.size(); ++v)
  {
    const double d2 = _spheres->_pt.distance_squared( c, vor_vertices[ cell0[v] ] );
    // don't trust a vertex that is more than 1/2 the neighbor distance away
    if (d2 > ignore_d2)
      continue;
    if ( d2 > max_d2 )
      max_d2 = d2;
  }
  // no nearby vertices?
  if (max_d2 == 0.)
    return 100.;
  
  const double beta = sqrt(max_d2) / _spheres->_pt.radius(c);
  return beta;
}


double MPS::estimate_beta( size_t sphere_index )
{
  // in high dimensions, only sampe the "sweet spot" where the theta directions are nicely spaced
  bool do_quadrant = num_dim() > 4;
  
  double beta = 0.; // small value
  // increment for theta. not sure what to pick
  // smaller is more accurate, but larger is faster
  // I think the true value is guaranteed to be at most
  //   max distance between starting point of rays * sqrt(d) / 2: about increment sqrt(d) / 2
  //  * max distance of found beta = 2 or so
  //  so PI * 0.03 means we're within PI*0.03 = 0.09 sqrt(d) / 2 about 0.1 * cos(30) = 0.05 or so in d=5
  // And likely to be much more
  //  const double dt = PI * 0.03; d=5,
  const double dt = PI * 0.05; // d=6
  //const double dt = PI * 0.025; // d=6
  // const double dt = PI * 0.001; // d=6
  const double *c = (*_spheres)[sphere_index];
  
  Trimming_Tool tt( num_dim() );

  double *g = tt.new_sphere();
  double *theta = tt.new_sphere();
  double *u = tt.new_sphere();
  double *p = tt.new_sphere();
  double *n = tt.new_sphere();
  
  const double A(0.5);

// size_t loop_count(0);
//  loop_count++;
//  if (loop_count == 27)
//    std::cout<<"debug me" << std::endl;
//  
  // for all polar directions theta
  if (do_quadrant)
    tt.first_quadrant_theta(theta);
  else
    tt.assign(theta, 0.);
  // debug
  // size_t theta_count(0);
  do
  {
    tt.theta_to_xyz(theta, u);
    double A_1(A), A_2(4.);
    
//    std::cout << "theta count " << theta_count << std::endl;
//    std::cout << theta_count << " ";
    
    // trim the spoke in the theta direction by the separating hyperspheres
    for (size_t si = 0; si < _spheres->num_real_spheres(); ++si)
    {
      if ( si == sphere_index )
        continue;
//      
//      if (loop_count == 27 && si == 24 && theta_count == 16)
//        std::cout<<"debug me" << std::endl;

      const double *s = (*_spheres)[si];
      bool is_s(true);
      tt.closest_ghost(c, s, g, is_s);
      
      tt.hyperplane_trim(c, u, s, A, A_1, A_2, p, n );
    
      // point at 0.5 can't be trimmed without violating the conflict sphere radius
      assert( !((A_2 < A * 0.99) || (A_1 > A * 1.01)) );
    
      // if we've already trimmed it shorter than the farthest-known beta,
      // then skip all other spheres and move on to the next angle
      if ( A_2 < beta )
        break;
    }
    
    // if we've trimmed by all the spheres and we're still farther than anything we've seen so far
    // then save the new beta
    if (A_2 > beta)
      beta = A_2;
    
//    theta_count++;
    
  } while ( do_quadrant ? tt.next_quadrant_theta(theta, dt) : tt.next_theta(theta, dt) );
  
  tt.delete_sphere(g);
  tt.delete_sphere(theta);
  tt.delete_sphere(u);
  tt.delete_sphere(p);
  tt.delete_sphere(n);

  return beta;
}
