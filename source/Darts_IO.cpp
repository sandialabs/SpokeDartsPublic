// Darts_IO.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Darts_IO.hpp"

#include <sstream>
#include <fstream>
#include <iomanip>
#include <functional>

#include "Spoke_Length.hpp"
#include "Ghost_Global_Container.hpp"
#include "Ghost_Spheres.hpp"
#include "Search_Structure.hpp"
#include "Sphere_Container_Array.hpp"

#include "Quality_Tool.hpp"

void  Darts_IO::ps_preamble(std::ostream &file)
{
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;
  
  // set font
  //  file << "/Times-Roman findfont\n" << "12 scalefont\n" << "setfont" << std::endl;
  file << "/Courier findfont\n" << "0.12 scalefont\n" << "setfont" << std::endl;

}

void Darts_IO::ps_seg(std::ostream &file, std::string &name, std::string &color, double linewidth)
{
  file << "/" << name << "      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " grestore" << std::endl;
  file << " " << color << " setrgbcolor" << std::endl;
  file << " " << linewidth << " setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
}

void Darts_IO::ps_filledcirc(std::ostream &file, std::string &name, std::string &fillcolor, std::string &linecolor, double linewidth)
{
  file << "/" << name <<"    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl; 
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " " << fillcolor << " setrgbcolor" << std::endl;
  file << " fill" << std::endl;
  file << " grestore" << std::endl;
  file << " " << linecolor << " setrgbcolor" << std::endl;
  file << " " << linewidth << " setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
}

void Darts_IO::ps_circ(std::ostream &file, std::string &name, std::string &linecolor, double linewidth)
{
  file << "/" << name << "    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " " << linecolor << " setrgbcolor" << std::endl;
  file << " " << linewidth <<" setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
}
  
void Darts_IO::ps_blackquad(std::ostream &file)
{
  file << "/blackquad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl; 
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " grestore" << std::endl;
  file << " 0 0 0 setrgbcolor" << std::endl;
  file << " 0.02 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
}

void Darts_IO::ps_quad_white(std::ostream &file)
{
  file << "/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl; 
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 1.0 setgray fill" << std::endl;
  file << " grestore" << std::endl;   
  file << "} def" << std::endl;
}

void Darts_IO::ps_quad_bold(std::ostream &file)
{
  file << "/quad_bold      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl; 
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.01 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
}


void Darts_IO::ps_primitives(std::ostream &file, const double base_r)
{
  // vocabulary:

  // see http://mcgraweng.com/Postscript%20Colors.pdf
  const int num_segs = 7;
  std::string  seg_names[num_segs] = {"redseg", "redsegthin", "greenseg", "bluesegfat", "blackseg", "bgsegfat", "orangeseg" };
  std::string seg_colors[num_segs] = { "1 0 0",      "1 0 0",    "0 1 0",      "0 0 1",    "0 0 0",    "0 1 1", "1 0.7 0" };
  double seg_widths[num_segs] = {        0.01,        0.005,       0.01,           1.0,       0.01,        1.0,    0.01 };
  for (size_t i = 0; i < num_segs; ++i)
    ps_seg( file, seg_names[i], seg_colors[i], seg_widths[i] * base_r);

  // f = filled circles, l = line around the circles
  const int num_fcircs = 6;
  std::string   fcirc_names[num_fcircs] = {"blackfcirc", "redfcirc", "bluefcirc", "greenfcirc", "greenflcirc", "bgreenflcirc" };
  std::string fcirc_fcolors[num_fcircs] = { "0 0 0",        "1 0 0",     "0 0 1",      "0 1 0",      "0 1 0",       "0 1 1" };
  std::string fcirc_lcolors[num_fcircs] = { "0 0 0",        "0 0 0",     "0 0 0",      "0 0 0",      "0 0 0",       "0 0 0" };
  double  fcirc_widths[num_fcircs] = {        0.00,           0.00,        0.00,         0.00,         0.01,          0.01 };
  for (size_t i = 0; i < num_fcircs; ++i)
    ps_filledcirc( file, fcirc_names[i], fcirc_fcolors[i], fcirc_lcolors[i], fcirc_widths[i] * base_r);

  const int num_circs = 2;
  std::string   circ_names[num_circs] = { "blacklcircfat", "blacklcirc"};
  std::string circ_lcolors[num_circs] = {         "0 0 0",     " 0 0 0"};
  double  circ_widths[num_circs]      = {           0.02,         0.01};
  for (size_t i = 0; i < num_circs; ++i)
    ps_circ( file, circ_names[i], circ_lcolors[i], circ_widths[i] * base_r);

  ps_blackquad(file);
  ps_quad_white(file);
  ps_quad_bold(file);
}


double Darts_IO::ps_translate_and_scale( std::ostream &file, Domain *domain )
{
  // region Retrieve bounding box, scale, and translate:
  double scale(1.);
  double shift_x(0.), shift_y(0.);
    
  const double Lx(domain->xmax()[0] - domain->xmin()[0]);
  const double Ly(domain->xmax()[1] - domain->xmin()[0]);        
    
  const double scale_x = 6.5 / Lx;
  const double scale_y = 9.0 / Ly;
    
  if (scale_x < scale_y) 
  {    
      scale = scale_x;                
      shift_x = 1.0 - domain->xmin()[0] * scale;
      shift_y = 0.5 * (11.0 - Ly * scale) - domain->xmin()[1] * scale;
  }
  else 
  {    
      scale = scale_y;        
      shift_x = 0.5 * (8.5 - Lx * scale) - domain->xmin()[0] * scale;
      shift_y = 1.0 - domain->xmin()[1] * scale;
  }    
  file << shift_x << " " << shift_y << " translate" << std::endl;
  return scale;
}

void Darts_IO::ps_draw_circs( circ_function *cf, std::ostream &file, Sphere_Container *sc)
{
  for (size_t i = sc->first(); i != sc->bad_sphere_index(); i = sc->next())
  {
    (*cf)( file, sc->_spheres[i] );
  }
}

void Darts_IO::ps_label_circs( std::ostream &file, double scale, Sphere_Container *sc)
{
  for (size_t i = sc->first(); i != sc->bad_sphere_index(); i = sc->next())
  {
    std::stringstream ss;
    ss << i;
    ps_draw_text(file, scale, ss.str(), sc->_spheres[i][0], sc->_spheres[i][1]);
  }
}

void Darts_IO::ps_draw_circle( std::ostream &file, const double scale, const double radius_factor, const std::string &circ_name, const size_t num_dim, const double* c )
{
  Point_Tool pt( num_dim );
  file << c[0] * scale << "  " << c[1] * scale << "  " << pt.radius(c) * scale * radius_factor << "  ";
  file << circ_name << std::endl;
}
void Darts_IO::ps_draw_dotted_circle( std::ostream &file, const double scale, const double radius_factor, const std::string &circ_name, const size_t num_dim, const double* c )
{
  ps_draw_circle( file, scale, radius_factor, circ_name, num_dim, c );
  ps_draw_circle( file, scale, radius_factor * 0.1, "blackfcirc", num_dim, c );
}

void Darts_IO::ps_draw_domain(  std::ostream &file, const double scale, Domain *domain )
{
  file << domain->xmin()[0] * scale << "  " << domain->xmin()[1] * scale << "  ";
  file << domain->xmax()[0] * scale << "  " << domain->xmin()[1] * scale << "  ";          
  file << domain->xmax()[0] * scale << "  " << domain->xmax()[1] * scale << "  ";          
  file << domain->xmin()[0] * scale << "  " << domain->xmax()[1] * scale << "  ";          
  file << "quad_bold"      << std::endl;
  
  if (draw_text())
  {
    std::stringstream ss_min, ss_max;
    ss_min << "(" << domain->xmin()[0] << "," << domain->xmin()[1] << ")";
    ss_max << "(" << domain->xmax()[0] << "," << domain->xmax()[1] << ")";
    ps_draw_text( file, scale, ss_min.str(), domain->xmin()[0], domain->xmin()[1] );
    ps_draw_text( file, scale, ss_max.str(), domain->xmax()[0], domain->xmax()[1] );
  }
  
  
}


void Darts_IO::ps_draw_text( std::ostream &file, const double scale, const std::string &text, double x, double y )
{
  const int offset(0);
  file << "newpath " << x * scale << " " << y * scale + offset << " moveto (" << text << ") show" << std::endl;
}


// I'm not sure what this is supposed to do
// if (false)
// {
//   // plot domain boundaries
//   const double DX = domain->xmax()[0] - domain->xmin()[0];
//   const double DY = domain->xmax()[1] - domain->xmin()[1];

//   file << (domain->xmin()[0] - DX) * scale << "  "  <<  domain->xmin()[1]        * scale << "  ";  
//   file << (domain->xmax()[0] + DX) * scale << "  "  <<  domain->xmin()[1]       * scale << "  ";          
//   file << (domain->xmax()[0] + DX) * scale << "  "  << (domain->xmin()[1] - DY) * scale << "  ";          
//   file << (domain->xmin()[0] - DX) * scale << "  "  << (domain->xmin()[1] - DY) * scale << "  ";          
//   file << "quad_white"      << std::endl;

//   file << (domain->xmin()[0] - DX) * scale << "  "  <<  domain->xmax()[1]        * scale << "  ";  
//   file << (domain->xmax()[0] + DX) * scale << "  "  <<  domain->xmax()[1]        * scale << "  ";          
//   file << (domain->xmax()[0] + DX) * scale << "  "  <<  (domain->xmax()[1] + DY) * scale << "  ";          
//   file << (domain->xmin()[0] - DX) * scale << "  "  <<  (domain->xmax()[1] + DY) * scale << "  ";          
//   file << "quad_white"      << std::endl;

//   file <<  domain->xmax()[0]       * scale << "  "  <<  (domain->xmin()[1] - DY) * scale << "  ";  
//   file << (domain->xmax()[0] + DX) * scale << "  "  <<  (domain->xmin()[1] - DY) * scale << "  ";          
//   file << (domain->xmax()[0] + DX) * scale << "  "  <<  (domain->xmax()[1] + DY) * scale << "  ";          
//   file <<  domain->xmax()[0]       * scale << "  "  <<  (domain->xmax()[1] + DY) * scale << "  ";          
//   file << "quad_white"      << std::endl;

//   file << (domain->xmin()[0] - DX) * scale << "  "  <<  (domain->xmin()[1] - DY) * scale << "  ";  
//   file <<  domain->xmin()[0]       * scale << "  "  <<  (domain->xmin()[1] - DY) * scale << "  ";          
//   file <<  domain->xmin()[0]       * scale << "  "  <<  (domain->xmax()[1] + DY) * scale << "  ";          
//   file << (domain->xmin()[0] - DX) * scale << "  "  <<  (domain->xmax()[1] + DY) * scale << "  ";          
//   file << "quad_white"      << std::endl;
// }

void Darts_IO::ps_draw_seg( std::ostream &file, const double scale, const std::string &segname, size_t num_dim, const double *c, double A_1, double A_2, const double *u )
{
  const double s0 = (c[0] + A_1 * u[0]) * scale;
  const double s1 = (c[1] + A_1 * u[1]) * scale;
  const double t0 = (c[0] + A_2 * u[0]) * scale;
  const double t1 = (c[1] + A_2 * u[1]) * scale;
  file << s0 << " " << s1  << " " << t0 << " " << t1 << " ";
  file << segname << std::endl;
}

void Darts_IO::ps_draw_plus( std::ostream &file, const double scale, double length, const std::string &segname, size_t num_dim, const double *p )
{
  const double c0 = p[0] * scale;
  const double c1 = p[1] * scale;
  const double r = length * scale;
  file << c0 - r << "  " << c1 << " " << c0 + r << "  " << c1 << " ";
  file << segname     << std::endl;
  file << c0 << "  " << c1 - r << " " << c0 << "  " << c1 + r << " ";
  file << segname     << std::endl;
}



void Darts_IO::plot_vertices_2d(std::string name, Spheres *spheres, Domain *domain, const double radius_factor,
                                const bool do_ghosts, const bool do_solid_disks, const bool do_rings_on_top,
                                const bool do_last_disk, const size_t active_sphere,
                                const double *dart,  const Spoke_Length *sl,  const double *u,  const double *p,
                                const double *dart2, const Spoke_Length *sl2, const double *u2, const double *p2)
{

  // open file
  // construct a name for it
  std::stringstream ss;
  ss << name << "_plot_2d";
  if (!do_ghosts)
    ss << "_no_ghosts";
  ss <<".ps";
  std::fstream file(ss.str().c_str(), std::ios::out);

  // preamble
  const double base_r = spheres->_pt.radius( (*spheres)[0] );
  ps_preamble(file);
  ps_primitives(file, base_r);
  const double scale = ps_translate_and_scale(file, domain);
  const double dot_r = 0.1 * base_r; // * radius_factor when used

  // get a global container, which only contains the ghosts if do_ghosts and the spheres are actually ghost_spheres
  Global_Container *gc = spheres && spheres->size() ?
    (do_ghosts ? Ghost_Global_Container::new_global_container(*spheres) : new Global_Container(*spheres))
    : 0;

  const bool is_active = (active_sphere != spheres->bad_sphere_index());

  // drawing
  // note later-drawn items may cover earlier-drawn items
  if (gc)
  {

    //== black filled circles, or green filled circles with black outline and center dot
    using namespace std::placeholders;

    circ_function cf;
    
    // draw solid filled circles, with or without lines around them
    if ( do_solid_disks )
      cf = std::bind( &Darts_IO::ps_draw_circle,          this, _1, scale, radius_factor,  "blackfcirc", gc->num_dim(), _2 );
    else
    {
      if (do_rings_on_top)
        cf = std::bind( &Darts_IO::ps_draw_circle,        this, _1, scale, radius_factor,  "greenfcirc", gc->num_dim(), _2 );
      else
        cf = std::bind( &Darts_IO::ps_draw_dotted_circle, this, _1, scale, radius_factor, "greenflcirc", gc->num_dim(), _2 );
    }
    ps_draw_circs( &cf, file, gc );

    // draw rings on top of all the filled disks
    if ( do_rings_on_top && !do_solid_disks )
    {
      cf = std::bind( &Darts_IO::ps_draw_circle,          this, _1, scale, radius_factor,  "blacklcirc", gc->num_dim(), _2 );
      ps_draw_circs( &cf, file, gc );
    }
    
    //== last inserted disk
    if (do_last_disk)
      ps_draw_dotted_circle( file, scale, radius_factor, "bgreenflcirc", spheres->num_dim(), (*spheres)[spheres->last_disk()] );

    //== active sphere
    if (is_active)
      ps_draw_dotted_circle( file, scale, radius_factor, "redfcirc", spheres->num_dim(), (*spheres)[active_sphere] );
  }
 
  // zzyk add drawing lines
  // zzyk add drawing points
  // draw points after segments, since the points may lie on top of them and are smaller
  
  //== draw dart
  bool is_first_dart  = is_active && ps_draw_dart( file, scale, dart, sl, u, p,     (*spheres)[active_sphere], dot_r, "redsegthin", "bluesegfat", "blackseg");
  bool is_second_dart = is_active && ps_draw_dart( file, scale, dart2, sl2, u2, p2,                         p, dot_r, "blackseg", "bgsegfat", "orangeseg" );

  //== plot domain boundaries
  ps_draw_domain( file, scale, domain);
  
  if (draw_text())
  {
    if (gc)
      ps_label_circs(file, scale, gc);
    
    // legend
    const double dy = -0.02;
    double y(dy);
    if (do_last_disk)
      ps_draw_text(file, scale, "last disk = cyan", 0, y);
    y += dy;
    if (is_active)
      ps_draw_text(file, scale, "front disk = red", 0, y);
    y += dy;
    if (is_first_dart)
    {
      ps_draw_text(file, scale, "first dart", 0, y);
      y += dy;
      ps_draw_text(file, scale, "  dart = blue line.  anchor = red plus.  picked point = black plus.", 0, y);
      y += dy;
    }
    if ( is_second_dart )
    {
      ps_draw_text(file, scale, "second dart", 0, y);
      y += dy;
      ps_draw_text(file, scale, "  dart = cyan line.  anchor = black plus.  picked point = red-green plus", 0, y);
      y += dy;
    }
  }

  file << "showpage" << std::endl;

  delete gc;
}

bool Darts_IO::ps_draw_dart(std::ostream &file, const double scale,
                  const double *dart,  const Spoke_Length *sl,  const double *u,  const double *p,
                  const double *dart_center, double dot_r,
                  const std::string &dart_seg, const std::string &sl_seg, const std::string &p_seg)
{
  
  // blue dot (small circle) at dart point
  const bool do_dart_line = sl && u;
  const bool do_dart_cross = dart;
  const bool do_p = p;
  if (do_dart_line)
  {
    ps_draw_seg( file, scale, sl_seg, 2, dart_center, sl->A_1, sl->A_2, u);
  }
  if (do_dart_cross)
  {
    const double cross_length = dot_r * 2.;
    ps_draw_plus(   file, scale, cross_length, dart_seg, 2, dart );
    ps_draw_circle( file, scale, dot_r,         "blackfcirc", 2, dart );
  }
  if (do_p)
  {
    const double cross_length = dot_r * 2.;
    ps_draw_plus(   file, scale, cross_length, p_seg, 2, p );
    ps_draw_circle( file, scale, dot_r,         "blackfcirc", 2, p );
  }
  return (do_dart_line || do_dart_cross || do_p);
}


void Darts_IO::RDF_histogram_data(std::ostream &out, Histogram &histogram)
{
  // bin width relative to sphere radius
  // #bins  #bins-for-one-radius  maxvalue averagevalue numpoints dimension
  out << histogram.bins.size() << " " <<  histogram.bins_per_rx << " " << histogram.maxbin << " " << histogram.bin_average << " " <<  histogram.num_spheres << " " << histogram.dim << " " << histogram.min_distance << std::endl;
  
  // relative number of elements in each bin
  for (unsigned int k = 0; k < histogram.bins.size(); k++)
  {
    out << std::setw(10);
    out << histogram.bins[k] << " ";
  }
  out << std::endl;

}


void Darts_IO::RDF_histogram(std::ostream &out, Histogram &histogram)
{
  double big_r = 0.1 * floor( (10. * histogram.bins.size()) / histogram.bins_per_rx );
  
  out << "Drawing histogram, r:" << histogram.rx;
  out << " distances between r and " << big_r << "r";
  
  if (histogram.look_periodic)
    out << ", periodic,";
  else
    out << ", window:[" << histogram.window_frame << "r, 1-" << histogram.window_frame << "r],";
  out << " dimension:" << histogram.dim << " num_spheres:" << histogram.num_spheres << std::endl;

  out << "Minimum encountered interdisk distance was " << histogram.min_distance << " r." << std::endl;
  
  const int maxstar = 80;
  const double f = maxstar / histogram.maxbin;

  // draw
  out << "num_source_points:" << histogram.num_windowed_points << std::endl;
  out << "num_distances:" << histogram.num_distances << " in histogram:" << histogram.num_histogram_distances << std::endl;
  for (unsigned int k = 0; k < histogram.bins.size(); k++)
  {
    // row header
    if ( k % histogram.bins_per_rx == 0)
      out << 1 + k/histogram.bins_per_rx << ":";
    else
      out << " :";
    if ( (k * 10) % histogram.bins_per_rx == 0)
      out << (k * 10 / histogram.bins_per_rx) - (k/histogram.bins_per_rx) * 10 << " ";
    else
      out << "  ";
    
    // floating value for later plotting
    out << std::setw(10);
    out << histogram.bins[k] * f;
    // characters
    for (unsigned int c = 0; c < histogram.bins[k] * f; c++)
    {
      if (c == floor(histogram.bin_average * f))
        out << "|";
      else
        out << "-";
    }
    out << std::endl;
  }
  
  out << "average bin height is " << histogram.bin_average *f << " out of " << maxstar << std::endl;
}

void Darts_IO::save_spheres_for_beta( Spheres *spheres, double ghost_frame, std::string name_prefix )
{
  // can't create all the ghosts for the center sphere if the frame is larger than 1/2
  if (ghost_frame > 0.49)
  {
    ghost_frame = 0.49;
    std::cout << "Warning, truncating frame at 0.49" << std::endl;
  }

  Ghost_Spheres *ghosts = spheres->cast_to_ghost();
  if ( !ghosts && ghost_frame>0.)
  {
    std::cerr << "Warning: save_spheres_for_beta called with non-ghosts and non-zero frame. Ghosts will not be written and beta may be innacurate." << std::endl;
  }
  if (!ghosts)
    save_spheres( spheres, Beta, name_prefix );
  else
  {
    double old_frame = ghosts->frame();
    ghosts->set_frame( ghost_frame ); 
    save_spheres( ghosts, Beta, name_prefix );
    ghosts->set_frame(old_frame);
  }
}

void Darts_IO::save_point(std::ostream &file, const double * p, size_t d, bool do_endl )
{
  assert( d );
  file << p[0];
  for (size_t i = 1; i < d; i++)
    file << " " << p[i];
  if (do_endl)
    file << std::endl;  
}
void Darts_IO::save_sphere(std::ostream &file, const double * s, size_t d, bool do_endl )
{
  // save_point(file, s, d, false);
  assert( d );
  for (size_t i = 0; i < d; i++)
    file << s[i] << " ";
  file << Point_Tool(d).radius(s);
  if (do_endl)
    file << std::endl;  
}


void Darts_IO::save_spheres( Spheres *spheres, Sphere_File_Format form, std::string name_prefix)
{

  size_t num_dim( spheres->num_dim() );
  const size_t num_spheres = spheres->num_real_spheres(); // ignores ghosts

  std::stringstream ss;
  bool do_radius(true); 
  bool do_bounding_box(true);
  bool do_ghost_frame(false);
  switch (form)
  {
    case PSA :
    {
      do_bounding_box = false;
      do_radius = false;
      ss  << name_prefix << "_psa_spheres" << ".txt";
    }
    break;
    case Beta :
    {
      do_bounding_box = false;
      do_radius = false;
      do_ghost_frame = true;
      ss  << name_prefix << "_beta_points" << ".txt";
    }
    break;

    case Spectrum :
    {
      do_bounding_box = false;
      do_radius = false;
      do_ghost_frame = false;
      ss << name_prefix << "_spectrum_points" << ".txt";
    }
    break;
      
    default:
    case Plain :
    {
      do_bounding_box=true;
      do_radius=true;
      ss << name_prefix << "_mps_spheres" << ".dat";
    }
    break;
  }

  // open file
  std::fstream file(ss.str().c_str(), std::ios::out);

  Global_Container gc(*spheres); // ignores any ghosts on purpose

  // count ghost points in the frame
  // caller must set frame
  size_t num_ghosts(0);
  Ghost_Spheres *ghosts = spheres->cast_to_ghost();
  if (do_ghost_frame)
  {
    if (ghosts)
    {
      // caller must set frame
      for (size_t i = gc.first(); i != gc.bad_sphere_index(); i = gc.next())
      {
        const double *s = (*spheres)[i];
        ghosts->make_periodic_sphere_copies( s );
        if (ghosts->periodic_copies_size())
          num_ghosts += ghosts->periodic_copies_size()-1;
      }
    }
  }

  // header
  switch (form)
  {
    case PSA:
    {
      // #spheres 
      file << num_spheres << std::endl;
      if (num_dim == 1)
      {
        std::cerr << "Warning: PSA file of 1-dimensional data not supported" << std::endl;
        file      << "Warning: PSA file of 1-dimensional data not supported" << std::endl;
        return;
      }
      if (num_dim != 2)
      {
        std::cerr << "Warning: PSA file contains only first two coordinates of " << num_dim << std::endl;
        // if we put this warning in the file, it will not be readable
        num_dim = 2;
      }
    }
    break;
    case Beta :
    {
      // dim
      // #spheres
      file << num_dim << std::endl << num_spheres + num_ghosts << std::endl;
    }
    break;

    case Spectrum :
    {
      ; // nada
    }
    break;
      
    default:
    case Plain:
    {      
      // dim #spheres
      file << num_dim << " " << num_spheres << std::endl;
    }
    break;
  }


  // precision for coordinates
  file.precision(15); file << std::fixed;

  // Bounding box
  if (do_bounding_box)
  {
    save_point( file, spheres->domain().xmin(), num_dim, true);
    save_point( file, spheres->domain().xmax(), num_dim, true);
  }

  // coordinates, radius
  if (do_radius)
  {
    for (size_t i = gc.first(); i != gc.bad_sphere_index(); i = gc.next())
      save_sphere( file, (*spheres)[i], num_dim, true);
  }
  else
  {
    for (size_t i = gc.first(); i != gc.bad_sphere_index(); i = gc.next())
      save_point( file, (*spheres)[i], num_dim, true);
  }

  // ghost points in the frame
  size_t num_written_ghosts(0);
  if (do_ghost_frame && ghosts && num_ghosts)
  {
    for (size_t i = gc.first(); i != gc.bad_sphere_index(); i = gc.next())
    {
      const double *s = (*spheres)[i];
      ghosts->make_periodic_sphere_copies( s );
      for ( size_t j = 0; j < ghosts->periodic_copies_size(); ++j )
      {
        if (!ghosts->is_real_index(j))
        {
          const double *g = ghosts->periodic_copies()[j];
          save_point( file, g, num_dim, true);
          ++num_written_ghosts;
        }
      }
    }
    assert( num_ghosts == num_written_ghosts );
  }

}

/*

void Darts_IO::load_spheres_for_beta( Spheres *spheres )
{
  std::fstream file("beta_points.txt", std::ios::in);
  
  // header;
  size_t num_dim(0), num_spheres(0);
  file << num_dim << std::endl << num_spheres + num_ghosts << std::endl;
  
  Domain *domain(0);
  
  
  Global_Container gc(*spheres); // ignores any ghosts on purpose
  
  // count ghost points in the frame
  // caller must set frame
  size_t num_ghosts(0);
  Ghost_Spheres *ghosts = spheres->cast_to_ghost();
  if (do_ghost_frame)
  {
    if (ghosts)
    {
      // caller must set frame
      for (size_t i = gc.first(); i != gc.bad_sphere_index(); i = gc.next())
      {
        const double *s = (*spheres)[i];
        ghosts->make_periodic_sphere_copies( s );
        if (ghosts->periodic_copies_size())
          num_ghosts += ghosts->periodic_copies_size()-1;
      }
    }
    else
      std::cerr << "Warning: saving file without ghosts, even though ghost_frame was non-zero." << std::endl;
  }
  

}
 */

