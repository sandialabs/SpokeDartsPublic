//  ScalingStudy.cpp
//  spokes
//
// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0
//
#include <sstream>
#include <fstream>
#include <iomanip>

#include "Scaling_Study.hpp"

#include "MPS_Spoke.hpp"
#include "MPS_Favored.hpp"
#include "MPS_Two_Spoke.hpp"
#include "MPS_Bridson.hpp"
#include "Timing_Stats.hpp"

// parameters
// standard
bool _all_searches_global(false);
bool _perm_ghosts(false);
// for very coarse domains with high dimensions, using array search only, do true, false

void write_header( std::ostream &file )
{
  file << "n" << ", " << "t" << ", " << "n_neighbor" << ", "  << "t_neighbor" << ", "  << "n_neighbor2" << ", "  << "t_neighbor2" << ", " << "r" << std::endl;
}

void write_n_t( std::ostream &file, MPS &mps )
{
  const size_t n = mps.output_size();
  const double t = mps.output_time();
  const double t_neighbor = mps.neighbor_time();
  const double n_neighbor = mps.neighbor_size();
  const double t_neighbor2 = mps.neighbor2_time();
  const double n_neighbor2 = mps.neighbor2_size();
  const double r = mps.output_radius();
  file << n << ", " << t << ", " << n_neighbor << ", "  << t_neighbor << ", "  << n_neighbor2 << ", "  << t_neighbor2 << ", " << r << std::endl;
}

double run_spokes(std::ostream &file, size_t dim, bool is_periodic, double r, Search_Factory::Search_Type search_type, size_t &output_n )
{
  MPS_Spoke mps;
  mps.create(dim, is_periodic, r, search_type, _all_searches_global, _perm_ghosts);
  write_n_t(file, mps);
  output_n = mps.output_size();
  return mps.output_time();
}

double run_favored(std::ostream &file, size_t dim, bool is_periodic, double r, Search_Factory::Search_Type search_type, size_t &output_n )
{
  MPS_Favored mps;
  mps.create(dim, is_periodic, r, search_type, _all_searches_global, _perm_ghosts);
  write_n_t(file, mps);
  output_n = mps.output_size();
  return mps.output_time();
}

double run_two(std::ostream &file, size_t dim, bool is_periodic, double r, Search_Factory::Search_Type search_type, size_t &output_n )
{
  MPS_Two_Spoke mps;
  mps.create(dim, is_periodic, r, search_type, _all_searches_global, _perm_ghosts);
  write_n_t(file, mps);
  output_n = mps.output_size();
  return mps.output_time();
}

double run_bridson(std::ostream &file, size_t dim, bool is_periodic, double r, Search_Factory::Search_Type search_type, size_t &output_n )
{
  MPS_Bridson mps;
  mps.create(dim, is_periodic, r, search_type, _all_searches_global, _perm_ghosts);
  write_n_t(file, mps);
  output_n = mps.output_size();
  return mps.output_time();
}

bool more_t( double t_spoke, double cutoff_t )
{
  // always false if cutoff_t is 0
  return t_spoke < cutoff_t;
  
}

bool more_n( size_t n_spoke, size_t cutoff_n )
{
  // always false if cutoff_n is 0
  return n_spoke < cutoff_n;
}


void scaling_by_n( std::string prefix, size_t dim, bool is_periodic,  Search_Factory::Search_Type search_type,
                  bool do_spoke = true, bool do_favored = false, bool do_two = false,
                  const double cutoff_t = 2., const size_t min_iter = 1, const size_t cutoff_n = 0,
                  bool do_bridson = false)
{
  // vary r, collect timings for each algorithm
  
  // Do at least min_iter runs of MPS   (set to 0 if no requirements)
  // Continue to do iterations until the time it takes is greater than cutoff_t              (set to 0 if no requirements)
  //   "  "  "                           number of samples produced is greater than cutoff_n (set to 0 if no requirements)
  //
  // Iterations will stop only when *all three* conditions are met
  
  Timing_Stats run_time;
  run_time.start_clock();

  // get file names
  std::stringstream file0ss, file1ss, file2ss, file3ss;
  
  // junk filenames in order to not overwrite any interesting data
  if (!do_spoke)
    file0ss << "00empty_";
  if (!do_favored)
    file1ss << "00empty_";
  if (!do_two)
    file2ss << "00empty_";
  if (!do_bridson)
    file3ss << "00empty_";
  
  file0ss << prefix << dim;
  file1ss << prefix << dim;
  file2ss << prefix << dim;
  file3ss << prefix << dim;
  if (is_periodic)
  {
    file0ss << "_P_";
    file1ss << "_P_";
    file2ss << "_P_";
    file3ss << "_P_";
  }
  else
  {
    file0ss << "_n_";
    file1ss << "_n_";
    file2ss << "_n_";
    file3ss << "_n_";
  }
  std::string search_s;
  switch (search_type)
  {
    case Search_Factory::ARRAY:
      search_s = "a";
      break;
    case Search_Factory::GRID:
      search_s = "g";
      break;
    case Search_Factory::RANGE:
      search_s = "r";
      break;
    case Search_Factory::TREE:
      search_s = "t";
      break;
    case Search_Factory::UNSPECIFIED:
      search_s = "u";
      break;
  }
  file0ss << search_s << "_0spoke.txt";
  file1ss << search_s << "_1favored.txt";
  file2ss << search_s << "_2two.txt";
  file3ss << search_s << "_3bridson.txt";
  
  // open files
  std::fstream file0( file0ss.str().c_str(), std::ios::out);
  std::fstream file1( file1ss.str().c_str(), std::ios::out);
  std::fstream file2( file2ss.str().c_str(), std::ios::out);
  std::fstream file3( file3ss.str().c_str(), std::ios::out);

  // header
  file0 << "0 " << dim << " " << is_periodic << " " << search_s <<   " spoke "     << dim << "d" << std::endl;
  file1 << "1 " << dim << " " << is_periodic << " " << search_s << " favored "     << dim << "d" << std::endl;
  file2 << "2 " << dim << " " << is_periodic << " " << search_s <<     " two "     << dim << "d" << std::endl;
  file3 << "3 " << dim << " " << is_periodic << " " << search_s <<     " bridson " << dim << "d" << std::endl;
  write_header(file0);
  write_header(file1);
  write_header(file2);
  write_header(file3);
  
  // number of seconds allowed
  // do one iteration past this value, which could take up to 4x as long if O(n^2)
//  const double cutoff = 540; // three minutes
//    const double cutoff = 180; // three minutes
//  const double cutoff = 24; // 2 seconds
//  const double cutoff = 4;
//  const double cutoff = 0.01; // testing
  
  // time taken last iter
  double t_spoke = 0.;
  double t_favored = 0.;
  double t_two = 0.;
  double t_bridson = 0.;
  
  // number of real samples produced last iter
  size_t n_spoke(0), n_favored(0), n_two(0), n_bridson(0);
  
  const double min_r0 = 1.; // MPS_Spoke::max_radius( is_periodic ); //  0.161;
  const double min_r1 = 1.; // MPS_Favored::max_radius( is_periodic ); // 0.102;
  const double min_r2 = 1.; // MPS_Two_Spoke::max_radius( is_periodic); // 0.071;
  const double min_r3 = 1.; // MPS_Bridson::max_radius( is_periodic); // 0.071;

//  double first_n = 10;
//  double r = 1. / ( pow( first_n, 1. / (double) (dim) ));
//  double r = 1. - (1. / (4. * (double) dim ) );
//  if ( r > min_r0 )
    double r = min_r0;
  
  size_t spoke_iter(0), favored_iter(0), two_iter(0), bridson_iter(0);
  
  for (int i=0; i<1000; ++i)
  {
    std::cout << std::endl << "scaling_by_n r = " << r  << std::endl;
    {
      if ( do_spoke && (r <= min_r0) && (more_t( t_spoke, cutoff_t) || more_n( n_spoke, cutoff_n ) || more_n( spoke_iter, min_iter)) )
      {
        t_spoke = run_spokes(file0, dim, is_periodic, r, search_type, n_spoke);
        spoke_iter++;
      }
      if ( do_favored && (r <= min_r1) && (more_t( t_favored, cutoff_t) || more_n( n_favored, cutoff_n ) || more_n( favored_iter, min_iter)))
      {
        t_favored = run_favored(file1, dim, is_periodic, r, search_type, n_favored);
        favored_iter++;
      }
      if ( do_two && (r <= min_r2) && (more_t( t_two, cutoff_t) || more_n( n_two, cutoff_n ) || more_n( two_iter, min_iter)) )
      {
        t_two = run_two(file2, dim, is_periodic, r, search_type, n_two);
        two_iter++;
      }
      if ( do_bridson && (r <= min_r3) && (more_t( t_bridson, cutoff_t) || more_n( n_bridson, cutoff_n ) || more_n( bridson_iter, min_iter)) )
      {
        t_bridson = run_bridson(file3, dim, is_periodic, r, search_type, n_bridson);
        bridson_iter++;
      }
    }
    
    // double n each time
    // this could quadruple time
    //    r *= 1./sqrt(2);
    r *= 1./ pow(2., 1./double(dim));
    
    if ( 1/r > 20000 )
      break;
  }
  
  file0.close();
  file1.close();
  file2.close();
  file3.close();  
  
  run_time.collect_stats();
  std::cout << "scaling_by_n runtime = " << run_time.cpu_time << std::endl;
}

void scaling_study_1()
{
  // spokes, favored, and two, for periodic domains, up to 7d
  // tree for speed

  size_t dim = 2;
  bool is_periodic = false;
  
  assert(0); // asserts should be compiled out for timing study
  
  bool do_spoke = true;
  bool do_favored = true;
  bool do_two = true;
  double cutoff = 540;
  size_t min_iter = 3;

  if (1)
  {
    for (dim = 2; dim <= 7; ++dim)
    {
      cutoff = 540;
      min_iter = 3;
      
      if (dim == 4)
        cutoff = 1400;
      if (dim == 5)
        cutoff = 2000;
      if (dim == 6)
      {
        min_iter = 2;
        cutoff = 3000;
      }
      if (dim == 7)
      {
        min_iter = 1;
        cutoff = 5000;
      }
//    scaling_by_n("scaling_study_1", dim, is_periodic, Search_Factory::ARRAY ); // array is brute force n^2, always worse than tree
//    scaling_by_n("scaling_study_1", dim, is_periodic, Search_Factory::GRID ); // some problems to debug for coarse
//    scaling_by_n("scaling_study_1", dim, is_periodic, Search_Factory::RANGE );
      is_periodic = false;
      scaling_by_n("scaling_study_1_", dim, is_periodic, Search_Factory::TREE,
                  do_spoke, false, false, cutoff, min_iter);
      is_periodic = true;
      scaling_by_n("scaling_study_1_", dim, is_periodic, Search_Factory::TREE,
                   do_spoke, do_favored, do_two, cutoff, min_iter);
    }
  }
}

void scaling_study_2()
{
  // periodic array spoke, favored, two
  bool is_periodic = true;
  size_t dim = 2;
  double cutoff = 1500;
  size_t min_iter = 3;
  for (dim = 2; dim <= 6; ++dim)
  {
    scaling_by_n("scaling_study_2_", dim, is_periodic, Search_Factory::ARRAY,
                 true, dim < 6, dim < 6, cutoff, min_iter);
  }
}

void scaling_study_3()
{
  // aperiodic array spoke, high dimensions
  bool is_periodic = false;
  size_t dim = 24;
  double cutoff = 2000;
  size_t min_iter = 1;
  for (dim = 24; dim > 1; --dim)
  {
    scaling_by_n("scaling_study_3_", dim, is_periodic, Search_Factory::ARRAY,
                 true, false, false, cutoff, min_iter);
  }
}

void scaling_study_6()
{
  // aperiodic array spoke, pushing stunt dimensions
  bool is_periodic = false;
  size_t dim = 25;
  double cutoff_t = 0.; // 1800;
  size_t min_iter = 1;
  size_t cutoff_n = 100000; // a hundred thousand
  for (size_t i = 1; i < 100; ++i )
  {
    scaling_by_n("scaling_study_6_", dim, is_periodic, Search_Factory::ARRAY,
                 true, false, false, cutoff_t, min_iter, cutoff_n );
    
    dim++;

    // double dimension each iteration
    // dim *= 2;
    
    // linear in dimension, so double the time cutoff each time
    // but the N in the k_2 d n N factor might be growing too fast...
    // cutoff *= 2.;
    
    // final n should be over 100,000 each time
  }
}


void scaling_study_4()
{
  // periodic spoke, pushing dimension 6
  bool is_periodic = true;
  double cutoff = 1500;

  size_t dim = 6;
  size_t min_iter = 4;
  scaling_by_n("scaling_study_4a_", dim, is_periodic, Search_Factory::TREE,
               false, true, true, cutoff, min_iter);
  scaling_by_n("scaling_study_4b_", dim, is_periodic, Search_Factory::ARRAY,
               true, true, true, cutoff, min_iter);

  min_iter = 3;
  dim = 7;
  scaling_by_n("scaling_study_4c_", dim, is_periodic, Search_Factory::ARRAY,
               true, true, true, cutoff, min_iter);
  

}

void scaling_study_5()
{
  // periodic tree spoke, pushing dimension
  bool is_periodic = true;
  size_t dim = 2;
  double cutoff = 1500;
  size_t min_iter = 3;
  for (dim = 7; dim <= 8; ++dim)
  {
    scaling_by_n("scaling_study_4a_", dim, is_periodic, Search_Factory::TREE,
                 true, false, false, cutoff, min_iter);
  }
  for (dim = 6; dim <= 8; ++dim)
  {
    scaling_by_n("scaling_study_4b_", dim, is_periodic, Search_Factory::TREE,
                 false, true, true, cutoff, min_iter);
  }
}

void scaling_study_7()
{
  // periodic tree spoke, pushing dimension
  bool is_periodic = true;
  size_t dim = 2;
  double cutoff = 1500;
  size_t min_iter = 4;
  for (dim = 6; dim <= 6; ++dim)
  {
    scaling_by_n("scaling_study_7_", dim, is_periodic, Search_Factory::TREE,
                 false, false, true, cutoff, min_iter);
  }
}


void scaling_study_8()
{
  // periodic spoke, pushing dimension 6
  bool is_periodic = true;
  double cutoff_t = 1;
  size_t min_iter = 1;
  size_t min_n = 200000;
  
  size_t dim = 7;
  for (dim = 7; dim < 30; ++dim)
  {
    scaling_by_n("scaling_study_8_", dim, is_periodic, Search_Factory::ARRAY,
                 true, false, false, cutoff_t, min_iter, min_n);
  }
}

void scaling_study_9()
{
  // periodic spoke, pushing dimension, using array search now
  bool is_periodic = true;
  double cutoff_t = 1;
  size_t min_iter = 1;
  size_t min_n = 32000; // 32000;
  
  _all_searches_global = true;
  _perm_ghosts = false;
  
  size_t dim = 2;
  for (dim = 2; dim < 9; ++dim)
  {
    scaling_by_n("scaling_study_9_", dim, is_periodic, Search_Factory::ARRAY,
                 true, true, true, cutoff_t, min_iter, min_n);
  }
}
void scaling_study_10()
{
  // periodic spoke, pushing dimension, using array search now
  bool is_periodic = true;
  double cutoff_t = 1;
  size_t min_iter = 1;
  size_t min_n = 32000; // 32000;
  
  _all_searches_global = true;
  _perm_ghosts = false;
  
  size_t dim = 2;
  for (dim = 9; dim < 13; ++dim)
  {
    scaling_by_n("scaling_study_10_", dim, is_periodic, Search_Factory::ARRAY,
                 true, true, true, cutoff_t, min_iter, min_n);
  }
}
void scaling_study_11()
{
  // periodic spoke, pushing dimension, using array search now
  bool is_periodic = true;
  double cutoff_t = 1;
  size_t min_iter = 1;
  size_t min_n = 32000; // 32000;
  
  _all_searches_global = true;
  _perm_ghosts = false;
  
  size_t dim = 2;
  for (dim = 13; dim < 16; ++dim)
  {
    scaling_by_n("scaling_study_11_", dim, is_periodic, Search_Factory::ARRAY,
                 true, true, true, cutoff_t, min_iter, min_n);
  }
}

void scaling_study_12()
{
  // periodic spoke, pushing dimension, using array search now
  bool is_periodic = true;
  double cutoff_t = 1;
  size_t min_iter = 1;
  size_t min_n = 32000; // 32000;
  
  _all_searches_global = true;
  _perm_ghosts = false;
  
  size_t dim = 2;
  for (dim = 16; dim < 64; dim += 2)
  {
    scaling_by_n("scaling_study_12_", dim, is_periodic, Search_Factory::ARRAY,
                 true, true, true, cutoff_t, min_iter, min_n);
  }
}


void scaling_study_bridson_1()
{
  bool is_periodic = true;
  double cutoff_t = 1;
  size_t min_iter = 1;
  size_t min_n = 32000; // 32000;
  
  _all_searches_global = true;
  _perm_ghosts = false;
  
  size_t dim = 2;
  for (dim = 2; dim < 7; ++dim)
  {
    min_n = 32000;
    scaling_by_n("scaling_study_bridson_1a_", dim, is_periodic, Search_Factory::ARRAY,
                 false, false, false, cutoff_t, min_iter, min_n, true);
    //    min_n = 200000;
    //    scaling_by_n("scaling_study_bridson_1t_", dim, is_periodic, Search_Factory::TREE,
    //                 false, false, false, cutoff_t, min_iter, min_n, true);
  }
  min_n = 32000;
  for (dim = 7; dim < 16; ++dim)
  {
    scaling_by_n("scaling_study_bridson_1b_", dim, is_periodic, Search_Factory::ARRAY,
                 false, false, false, cutoff_t, min_iter, min_n, true);
  }
  for (dim = 16; dim < 32; dim += 2)
  {
    scaling_by_n("scaling_study_bridson_1c_", dim, is_periodic, Search_Factory::ARRAY,
                 false, false, false, cutoff_t, min_iter, min_n, true);
  }
  
}

void scaling_study_bridson_2()
{
  bool is_periodic = true;
  double cutoff_t = 1;
  size_t min_iter = 1;
  size_t min_n = 32000; // 32000;
  
  _all_searches_global = true;
  _perm_ghosts = false;
  
  size_t dim = 2;
  for (dim = 2; dim < 7; ++dim)
  {
    min_n = 32000;
    scaling_by_n("scaling_study_bridson_1a_", dim, is_periodic, Search_Factory::ARRAY,
                 false, false, false, cutoff_t, min_iter, min_n, true);
    //    min_n = 200000;
    //    scaling_by_n("scaling_study_bridson_1t_", dim, is_periodic, Search_Factory::TREE,
    //                 false, false, false, cutoff_t, min_iter, min_n, true);
  }
  min_n = 32000;
  for (dim = 7; dim < 16; ++dim)
  {
    scaling_by_n("scaling_study_bridson_1b_", dim, is_periodic, Search_Factory::ARRAY,
                 false, false, false, cutoff_t, min_iter, min_n, true);
  }
  for (dim = 16; dim < 32; dim += 2)
  {
    scaling_by_n("scaling_study_bridson_1c_", dim, is_periodic, Search_Factory::ARRAY,
                 false, false, false, cutoff_t, min_iter, min_n, true);
  }
  
}


void beta_study( size_t num_dim, bool estimate_beta, size_t layer_limit )
{
  const double r = 0.01;
  const bool is_periodic = false;
  
  size_t num_trials(100000); // 100000;
  
  //  num_trials = 10000; // 10,000 for d=4 takes 2.3 days
  // num_trials = 100;
  // num_trials -= 58313;
  
  //  num_trials = 10000; // 10,000 for d=5 takes 8 days with exact beta
  // num_trials = 10000; // 10,000 for d=6 takes 3 days with estimated beta
  
  //  num_trials = 10; // 10 quick test
//  num_trials = 1; // 10 quick test
  
  //  Search_Factory::Search_Type search_type = Search_Factory::ARRAY;
  //  bool all_searches_global = true;
  //  bool perm_ghosts = false;
  
  Search_Factory::Search_Type search_type = Search_Factory::TREE;
  bool all_searches_global = false;
  bool perm_ghosts = false;
  
  const bool do_spoke =
    false;
  // true;
  const bool do_two =
    // false;
    true;
  const bool do_bridson =
    false;
    // true;
  const bool do_favored =
    false;
    // true;
  
  //    Random::random_instance.initiate_random_generator(230950923095203985);
  // Random::random_instance.initiate_random_generator(345093405903495);
  //  Random::random_instance.initiate_random_generator(49360349817350423);
  //  Random::random_instance.initiate_random_generator(5406932725678023456);
  //    Random::random_instance.initiate_random_generator(345346905409827103);
  //    Random::random_instance.initiate_random_generator(9327678023456233456);
  
  
  // spoke
  if (do_spoke)
  {
    MPS_Spoke mps;
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    mps._layer_limit = layer_limit;
    // turn off all debugging and other file output in MPS_Spoke.cpp
    
    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      //    mps_spoke.report_central_beta("spoke");
      mps.report_central_beta("spoke", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with spoke" << std::endl;
  }
  
  // two
  if (do_two)
  {
    MPS_Two_Spoke mps;
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    mps._layer_limit = layer_limit;

    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      mps.report_central_beta("two", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with two" << std::endl;
  }
  
  // bridson
  if (do_bridson)
  {
    MPS_Bridson mps;
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    mps._layer_limit = layer_limit;
    
    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      mps.report_central_beta("bridson", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with bridson" << std::endl;
  }
  
  
  // favored - this is the slowest
  if (do_favored)
  {
    MPS_Favored mps;
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    mps._layer_limit = layer_limit;
    // turn off all debugging and other file output in MPS_favored.cpp
    
    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      mps.report_central_beta("favored", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with favored" << std::endl;
  }
} // beta_study



void beta_study2( size_t num_dim, bool estimate_beta )
{
  const double r_factor = 1.5; // r_factor == 1 gives 32000 points
  const bool is_periodic = true;
  
  size_t num_trials(100); // 100000;
  
  //  num_trials = 10000; // 10,000 for d=4 takes 2.3 days
  // num_trials = 100;
  // num_trials -= 58313;
  
  //  num_trials = 10000; // 10,000 for d=5 takes 8 days with exact beta
  num_trials = 5000;
  
  //  num_trials = 10; // 10 quick test
  // num_trials = 1; // 10 quick test
  
  Search_Factory::Search_Type search_type = Search_Factory::ARRAY;
  bool all_searches_global = true;
  bool perm_ghosts = false;
  
  //  Search_Factory::Search_Type search_type = Search_Factory::TREE;
  //  bool all_searches_global = false;
  //  bool perm_ghosts = false;
  
  const bool do_spoke =
   false;
  // true;
  const bool do_two =
   false;
  //true;
  const bool do_bridson =
    false;
  // true;
  const bool do_favored =
   // false;
   true;
  
  //    Random::random_instance.initiate_random_generator(230950923095203985);
  // Random::random_instance.initiate_random_generator(345093405903495);
  //  Random::random_instance.initiate_random_generator(49360349817350423);
  //  Random::random_instance.initiate_random_generator(5406932725678023456);
  //    Random::random_instance.initiate_random_generator(345346905409827103);
  //    Random::random_instance.initiate_random_generator(9327678023456233456);
  
  
  // spoke
  if (do_spoke)
  {
    MPS_Spoke mps;
    
    const double r = 0.17 * r_factor;
    //0.17 for  d = 6 -> 32,000
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    // turn off all debugging and other file output in MPS_Spoke.cpp
    
    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      //    mps_spoke.report_central_beta("spoke");
      mps.report_central_beta("spoke", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with spoke" << std::endl;
  }
  
  // two
  if (do_two)
  {
    MPS_Two_Spoke mps;
    
    const double r = 0.09 * r_factor; //0.09 for  d = 6 -> 32,000
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    
    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      mps.report_central_beta("two", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with two" << std::endl;
  }
  
  // bridson
  if (do_bridson)
  {
    MPS_Bridson mps;
    
    const double r = 0.17 * r_factor;
    //0.17 for  d = 6 -> 32,000
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    
    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      mps.report_central_beta("bridson", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with bridson" << std::endl;
  }
  
  
  // favored - this is the slowest
  if (do_favored)
  {
    MPS_Favored mps;
    
    const double r = 0.14 * r_factor; //0.14 for  d = 6 -> 32,000
    
    mps._do_central_beta = true;
    mps._do_plain = false;
    mps._do_beta = false;
    mps._do_performance = false;
    // turn off all debugging and other file output in MPS_favored.cpp
    
    // empt beta.txt file
    //  /* int ok = */ system( "rm beta.txt; touch beta.txt" );
    
    for (size_t i = 0; i < num_trials; ++i)
    {
      // report progress
      size_t j = i % (1 + num_trials / 10);
      if ( j == 0 )
      {
        j = (i*10) / num_trials;
        std::cout << j <<  " ";
      }
      
      // random number generator should continue where it left off, not reseed
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
      mps.report_central_beta("favored", estimate_beta);
      mps.clear();
    }
    std::cout << std::endl << "done with favored" << std::endl;
  }
}




