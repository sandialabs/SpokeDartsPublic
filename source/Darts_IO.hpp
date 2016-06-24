// Darts_IO.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef DARTS_IO_HPP
#define DARTS_IO_HPP

#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <functional>
class Domain;
class Spoke_Length;
class Timing_Stats;
#include "Spheres.hpp"
class Sphere_Container;
class Search_Structure;
class Histogram;

//typedef std::function<void (const double *)> circ_function;
typedef std::function<void (std::ostream &file, const double *)> circ_function;

class Darts_IO
{
public:

  // ===========================
  // idiom for files / std::cout
  //============================
  //
  // std::stringstream name;
  // name << "mps_spheres_" << total_number_of_flats << "_RDF.txt";
  // ostream &out = new_ostream( use_stdout, name );
  // delete_ostream( out );
  //
  // ostream &out = new_ostream( true );
  // delete_ostream( true, out ); // optional
  //
  std::ostream &new_ostream( bool use_stdout )
  {
    return new_ostream_s( use_stdout, std::string() );
  }
  std::ostream &new_ostream( bool use_stdout, std::stringstream &ss )
  {
    return new_ostream_s( use_stdout, ss.str() );
  }
  std::ostream &new_ostream_s( bool use_stdout, std::string s = std::string() )
  {
    if (use_stdout || s.empty())
      return std::cout;
    return new_file_ostream(s);
  }
  std::ostream &new_file_ostream( std::stringstream &ss )
  {
    return new_file_ostream(ss.str());
  }
  std::ostream &new_file_ostream( std::string s )
  {
    std::ofstream *file_out = new std::ofstream( s.c_str(), std::ios::out );
    return *file_out;
  }
  void delete_ostream( bool use_stdout, std::ostream &f = std::cout )
  {
    if (use_stdout || ( &f == &std::cout ) )
      return;
    delete &f;
  }
  
  std::string loop_string( size_t loop_count, std::string suffix = std::string() )
  {
    std::stringstream ss;
    ss << loop_count;
    if ( suffix.size() )
      ss << "_" << suffix;
    return ss.str();
  }


  //==========================
  // settings
  //==========================
  void set_draw_text(bool new_val) {_draw_text = new_val;}
  const bool draw_text() {return _draw_text;}
  
  Darts_IO() : _draw_text(false) {;}
  
  //==========================
  // data i/o
  //==========================
  enum Sphere_File_Format {Plain, PSA, Beta, Spectrum};
  void save_spheres( Spheres *spheres, 
    Sphere_File_Format form = Darts_IO::Plain, std::string name_prefix = std::string() );
  void save_spheres_for_psa( Spheres *spheres, std::string name_prefix = std::string())
  { save_spheres( spheres, Darts_IO::PSA, name_prefix ); }
  // For anjul's matlab code
  void save_spheres_for_spectrum( Spheres *spheres, std::string name_prefix = std::string())
  { save_spheres( spheres, Darts_IO::Spectrum, name_prefix ); }
  void save_spheres_plain(Spheres *spheres, std::string name_prefix = std::string() )
  { save_spheres( spheres, Darts_IO::Plain, name_prefix ); }
  void save_spheres_for_beta( Spheres *spheres, double ghost_frame, std::string name_prefix = std::string() );
  
  //==========================
  // ps plots
  //==========================
  void plot_vertices_2d_allblack(std::string name, Spheres *spheres, Domain *domain, const double radius_factor)
  {
    const bool do_ghosts = false;
    plot_vertices_2d( name, spheres, domain, radius_factor, do_ghosts, true, false, false,
                     spheres->bad_sphere_index(), 0, 0, 0, 0, 0, 0, 0, 0 );
  }
  
  void plot_vertices_2d_rings(std::string name, Spheres *spheres, Domain *domain, const double radius_factor,
                               const size_t active_sphere = Spheres::bad_sphere_index() )
  {
    const bool do_ghosts = true;
    plot_vertices_2d( name, spheres, domain, radius_factor, do_ghosts, false, true,
                      true, active_sphere, 0, 0, 0, 0, 0, 0, 0, 0 );
  }
  
  void plot_vertices_2d_dart(std::string name, Spheres *spheres, Domain *domain, const double radius_factor,
                             const size_t active_sphere, const double *dart, const Spoke_Length *sl, const double *u, const double *p = 0)
  {
    const bool do_ghosts = true;
    plot_vertices_2d( name, spheres, domain, radius_factor, do_ghosts, false, false,
                     true, active_sphere, dart, sl, u, p, 0, 0, 0, 0);
  }
  void plot_vertices_2d_dart2(std::string name, Spheres *spheres, Domain *domain, const double radius_factor,
                              const size_t active_sphere,
                              const double *dart,  const Spoke_Length *sl,  const double *u,  const double *p,
                              const double *dart2, const Spoke_Length *sl2, const double *u2, const double *p2)
  {
    const bool do_ghosts = true;
    plot_vertices_2d( name, spheres, domain, radius_factor, do_ghosts, false, false,
                     true, active_sphere, dart, sl, u, p, dart2, sl2, u2, p2);
  }
  
  // include the last parameter if you want the active disk in red
  void plot_vertices_2d_nodart(std::string name, Spheres *spheres, Domain *domain, const double radius_factor,
                               const size_t active_sphere = Spheres::bad_sphere_index() )
  {
    const bool do_ghosts = true;
    plot_vertices_2d( name, spheres, domain, radius_factor, do_ghosts, false, false,
                      true, active_sphere, 0, 0, 0, 0, 0, 0, 0, 0 );
  }

 
  void plot_vertices_2d(std::string name, Spheres *spheres, Domain *domain, const double radius_factor,
                        const bool do_ghosts, const bool do_solid_disks, const bool do_rings_on_top,
                        const bool do_last_disk, const size_t active_sphere,
                        const double *dart,  const Spoke_Length *sl,  const double *u,  const double *p,
                        const double *dart2, const Spoke_Length *sl2, const double *u2, const double *p2);
//      const size_t num_segments, double* st, double* end, double* line_p, 
//       bool scale_down, size_t iflat,  
//                                     std::vector< double *> *points)

// old
  // void save_spheres(size_t num_dim, size_t num_spheres, double** spheres, double* xmin, double* xmax, size_t number_flats, double cpu_time);
  
  // void save_spheres_for_beta(size_t num_dim, size_t num_spheres, double** spheres);
  // void save_spheres_for_psa(size_t num_dim, size_t num_spheres, double** spheres, double* xmin, double* xmax, size_t number_flats, double cpu_time,
  //   const size_t num_ghosts = 0, const size_t *ghost_sphere_indices = NULL);
  // void report_timing(size_t num_spheres,size_t num_flats, size_t num_successive_misses,double cpu_time);

  // plot_vertices_2d(mylibs._num_inserted_spheres, mylibs._spheres, 1, st, end, line_p, xmin, xmax, scd.iactive, false, 0, scd.dart,
  //                        0, 0, 0, &points); //zzyk
  // void Balloon_Darts::plot_vertices_2d(const size_t num_spheres, double** spheres, const size_t num_segments, double* st, double* end, double* line_p, const double* xmin, const double *xmax,
  //                                  const size_t active_sphere, bool scale_down, size_t iflat, const double* dart, const bool do_black_half_disks,
  //                                    const size_t num_ghosts, const size_t *ghost_sphere_indices,
  //                                    std::vector< double *> *points);


  // ===========
  // text plots
  // ===========
  // write to file RDF_histogram_id_name.txt  human readable
  void RDF_histogram(std::string name, Histogram &histogram)
  {
    std::stringstream fname;
    fname << name << "_RDF_histogram.txt";
    std::ofstream file_out( fname.str().c_str(), std::ios::out );
    RDF_histogram(file_out, histogram);
  }
  // human readable output
  void RDF_histogram(std::ostream &out, Histogram &histogram);
  // hint:
  // Global_Search search(spheres, KD_Tree or Sphere_Array);
  // search = Ghost_Global_Search::new_global_search(spheres, KD_Tree or Sphere_Array);

  // values for further processing
  void RDF_histogram_data(std::string name, Histogram &histogram)
  {
    std::stringstream fname;
    fname << name << "_RDF_histogram_data.txt";
    std::ofstream file_out( fname.str().c_str(), std::ios::out );
    RDF_histogram_data(file_out, histogram);
  }
  void RDF_histogram_data(std::ostream &out, Histogram &histogram);

protected:
  
  // settings
  bool _draw_text;
  
  // write coordinates (and radius)
  void save_point (std::ostream &file, const double * p, size_t d, bool do_endl = true );
  void save_sphere(std::ostream &file, const double * s, size_t d, bool do_endl = true );

  // postscript primitives
  void ps_preamble(std::ostream &file);
  void ps_primitives(std::ostream &file, const double base_r = 1.);
  void ps_seg( std::ostream &file, std::string &name, std::string &color, double linewidth);
  void ps_filledcirc(std::ostream &file, std::string &name, std::string &fillcolor, std::string &linecolor, double linewidth);
  void ps_circ(std::ostream &file, std::string &name, std::string &linecolor, double linewidth);
  void ps_blackquad(std::ostream &file);
  void ps_quad_white(std::ostream &file);
  void ps_quad_bold(std::ostream &file);
  // returns scale
  double ps_translate_and_scale( std::ostream &file, Domain *domain );
  
  // some drawing function with everythign but the particular disk c filled in
  // e.g.
  //     circ_function cf = std::bind( &Darts_IO::ps_draw_circle, this, _1, scale, radius_factor, "blackfcirc", gc->num_dim(), _2 ) :
  void ps_draw_circs( circ_function *cf, std::ostream &file, Sphere_Container *sc);
  void ps_label_circs( std::ostream &file, double scale, Sphere_Container *sc);


  void ps_draw_circle(  std::ostream &file, const double scale, const double radius_factor, const std::string &circ_name, const size_t num_dim, const double* c);
  void ps_draw_dotted_circle(  std::ostream &file, const double scale, const double radius_factor, const std::string &circ_name, const size_t num_dim, const double* c);
  void ps_draw_domain(  std::ostream &file, const double scale, Domain *domain );
  void ps_draw_plus( std::ostream &file, const double scale, double length, const std::string &segname, size_t num_dim, const double *p );
  void ps_draw_seg( std::ostream &file, const double scale, const std::string &segname, size_t num_dim, const double *c, double A_1, double A_2, const double *u );
  void ps_draw_text( std::ostream &file, const double scale, const std::string &text, double x, double y );

  bool ps_draw_dart(std::ostream &file, const double scale,
                    const double *dart,  const Spoke_Length *sl,  const double *u,  const double *p,
                    const double *dart_center, // e.g. sphere center
                    double dot_r, // draw size radius of a point
                    const std::string &dart_seg, const std::string &sl_seg, const std::string &p_seg);

};

#endif