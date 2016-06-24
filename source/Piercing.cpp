// Piercing.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Piercing.hpp"
#include "Sphere_Container.hpp"

bool Piercing_Tool::two_piercings(const double A, const double B, const double x, const double y,
                                 Piercings &piercings, size_t &depth_at_A)
{
  assert( x <= y );

  // sphere covers A ?
  if ( x < A && y > A )
  {
    depth_at_A++;
    if ( y > B)
      return true;
  }
 
  // x between A and B ?
  if ( x >= A && x <= B )
    piercings.push_back(Piercing(x, true));

  // y between A and B ?
  if ( y >= A && y <= B )
    piercings.push_back(Piercing(y, false));  

  // whole segment not covered
  return false;
}

bool Piercing_Tool::uncovered_segments(Piercings &piercings, const size_t &depth_at_A, double A, double B,
  std::vector<double> &segments, double &sum_lengths, double min_length )
{
  sort_piercings( piercings );

  // traverse to find uncovered segment lengths
  size_t current_depth = depth_at_A;
  sum_lengths = 0.;
  segments.clear();
  double segA = A; // start of last uncovered segment
  for ( size_t i = 0; i < piercings.size(); ++i )
  {
    if (piercings[i].forward)
    {
      ++current_depth;
      // did I just finish an interval?
      if ( current_depth == 1)
      {
        assert( piercings[i].p - segA >= 0.);
        const double segB = piercings[i].p;
        double length = segB - segA;
        if (length > min_length)
        {
          sum_lengths += length;
          segments.push_back( segA );
          segments.push_back( segB );
        }
      }
    }
    else
    {
      assert( current_depth > 0 );
      --current_depth;
      // am I starting a new interval?
      if (current_depth == 0)
      {
        segA = piercings[i].p;
      }
    }
  }
  // is last segment open ended?
  if (current_depth == 0)
  {
    double length = B - segA;
    if (length > min_length)
    {
      sum_lengths += length;
      segments.push_back( segA );
      segments.push_back( B ); 
    }
  }
  return !segments.empty();
}


bool Piercings_Tool::gather_piercings( Sphere_Container *pierced_spheres,
    const double *p, const double *u, const double A, const double B, 
    Piercings *piercings_p, size_t *depth_at_A_p,
    Piercings *piercings_m, size_t *depth_at_A_m,
    bool quit_if_segment_touched )
{
  //for all spheres, 
  //  find the piercings, 
  //  add x,y to piercings_p and and -y,-x to piercings_m,
  double x, y;
  bool covered_p = false;
  bool covered_m = false;
  for (size_t i = pierced_spheres->first(); i != pierced_spheres->bad_sphere_index(); i = pierced_spheres->next())
  {
    const double *s = _spheres[i];
    _tt.points_of_line_piercing_sphere( p, u, s, x, y );
    covered_p = covered_p ||
      two_piercings( A,B, x,y, *piercings_p, *depth_at_A_p);
    if (piercings_m)
    {
      covered_m = covered_m ||
        two_piercings( A,B, -y,-x, *piercings_m, *depth_at_A_m);
    }
    // quit? 
    if ( (covered_p && covered_m) ||
      ( quit_if_segment_touched && piercings_p->size() && (piercings_m == 0 || piercings_m->size())) )
      return true;
  }
  return false;
}

bool Piercings_Tool::pierce_domain( 
    const double *xmax, const double *xmin,
    const double *p, const double *u, const double A, const double B, 
    Piercings *piercings_p, size_t *depth_at_A_p,
    Piercings *piercings_m, size_t *depth_at_A_m,
    bool quit_if_segment_touched)
{
  bool touched_p(true), covered_p(true);
  pierce_domain_oneside( xmax, xmin, p, u, A, B, piercings_p, depth_at_A_p, touched_p, covered_p );
  bool touched_m(true), covered_m(true);
  if ( depth_at_A_m )
  {
    pierce_domain_oneside( xmax, xmin, p, u, -A, -B, piercings_m, depth_at_A_m, touched_m, covered_m );
  }
  return ( covered_m && covered_p ) || (quit_if_segment_touched && (touched_p || touched_m));
}

void Piercings_Tool::pierce_domain_oneside( 
    const double *xmax, const double *xmin,
    const double *p, const double *u, const double A, const double B, 
    Piercings *piercings, size_t *depth_at_A,
    bool &touched, bool &covered)
{
  touched = false;
  covered = false; 
  double t(B);
  _tt.trim_by_domain( p, u, t, xmax, xmin );
  if (t <= A)
  {
    // whole positive interval is covered
    piercings->clear();
    ++(*depth_at_A);
    covered = true; 
    touched = true;
  }
  else if (t < B)
  {
    piercings->push_back(Piercing(t, true));
    touched = true;
  }
}
    
