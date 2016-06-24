// Spoke_Darts.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Spoke_Darts.hpp"
#include "Crossing_Tool.hpp"
#include "Sphere_Container.hpp"
#include "Timing_Stats.hpp" // Wheel_Stats
#include "Random.hpp"
#include "Search_Structure.hpp"


bool Spoke_Darts::wheel(const double *c, const double rc, double * u, double * dart,
           Sphere_Container *trim_spheres, const double radius_factor)
{
  //, double * xmin, double * xmax )
  Spoke_Darts_Workspace *ws = &spoke_darts_workspace;
  double *v = ws->v;

  _rng->random_orthonormal_vector(u, v);

  double theta_start, theta_end;
  if ( pick_segment_crossing( c, rc, u, v, trim_spheres, radius_factor, theta_start, theta_end, 0.))
  { 
    // pick a point of the uncovered interval
    double theta = _rng->random_middle(theta_start, theta_end);
    
    double a, b;
    theta_to_xy(theta, a,b);
    
    // new u
    axpby( ws->w, a, u, b, v );
    assign(u, ws->w);

    // new dart direction and magnitude, from c
    axpy( dart, rc, u, c );

    return true;
  }
  // depth zero never reached, whole circle was covered
  return false;
}

bool Spoke_Darts::pick_segment_crossing(const double *c, const double rc, double *u, double *v,
                                        Sphere_Container *trim_spheres, const double radius_factor,
      double &theta_start, double &theta_end, const double min_length )
{
  // outline
  // generate a plane, by picking a second dart.
  // generate a circle, parameterize it by t from 0 to 2pi
  // find the stop and start points of segments covered by ext_neighbor disks, by t, and the depth of the zero point.
  // sort the list of points
  // traverse the list, adding and subtracting to the depth,
  // for the longest uncovered arc, choose a new anchor point, favoring the middle of the arc
  // return the anchor point in dart
  
  // define the circle: the plane spanned by u and v
  // unit direction vectors from sphere center

  // collect the crossings
  size_t depth (0);
  Spoke_Darts_Workspace *ws = &spoke_darts_workspace;
  Crossings &crossings = ws->crossings;
  Crossings_Tool ct(_spheres);
  ct.gather_crossings(c, rc, u, v, trim_spheres, radius_factor,
      crossings, depth);

  assert(_domain.is_periodic());
  // clipping circles by the domain is not implemented yet. 
  // It is a bit of math and we don't really need it for the paper
  //   http://math.stackexchange.com/questions/213545/solving-a-trignometric-equation-of-form-a-sin-x-b-cos-x-c
  //  if (! _is_periodic )
  //    gather_crossings_axis( ci, A, u, v, xmin, xmax, crossings, depth);

  double sum_lengths(0);
  std::vector<double> &arcs = ws->arcs;
  arcs.clear();

  ct.uncovered_arcs(crossings, depth, arcs, sum_lengths, min_length);
  if ( arcs.empty() )
    return false;

  //zzyk to do: this last choice of how to pick an arc could be changed to something else
  // the above could be moved to one function, and the following to another
  ct.longest_interval( arcs, theta_start, theta_end);
  return true;

}


bool Spoke_Darts::pick_segment_piercing(const double* p, double* u, const double A, const double B, 
                                        Sphere_Container *trim_spheres,
                                        double & seg_start, double & seg_end,
                                        bool full_segment_only, double min_length)
{
  // outline
  // create a list of all crossings, where a crossing is where a sphere intersects the ray u from ci
  //   return false early if doing full_segment_only selection and there is one crossing
  // find the stop and start points of segments covered by ext_neighbor disks, by t, and the depth of the original dart point.
  //   sort the list of crossing points
  //   traverse the list, adding and subtracting to the depth
  // return false if no uncovered segments
  // pick one of the uncovered segments, uniform by segment length
  // return that segment as t_start, t_end
  // the caller will pick a point depending on the current rule, whether favored or unfavored
  // return true, an uncovered segment was found
  
  // re-used workspace for no dynamic memory
  Spoke_Darts_Workspace *ws = &spoke_darts_workspace;
  Piercings &piercings_p( ws->piercings_p );
  Piercings &piercings_m( ws->piercings_m );

  piercings_p.clear();  piercings_m.clear();
  size_t depth_p(0), depth_m(0); // depth of A and -A inside spheres
 

  Piercings_Tool pt(_spheres);
  bool quit = pt.gather_piercings( trim_spheres,
    p, u, A, B, 
    & piercings_p, & depth_p,
    & piercings_m, & depth_m,
    full_segment_only);

  if (quit)
    return false;

  // add dummy piercing coverage for non-periodic domain boundaries
  if (_domain.is_periodic())
  {
    quit = pt.pierce_domain( _domain.xmax(), _domain.xmin(), p, u, A, B,
        & piercings_p, & depth_p,
        & piercings_m, & depth_m,
        full_segment_only);
    if (quit)
      return false;
  }
  
  if (full_segment_only)
  {
    seg_start = A;
    seg_end = B;
    if (piercings_p.empty() && depth_p == 0)
    {
      return true;

    }
    else if (piercings_m.empty() && depth_m == 0)
    {
      multiply(u, -1.);
      return true;
    }
    return false;
  }

  std::vector<double> &segments_p( ws->segments_p );
  std::vector<double> &segments_m( ws->segments_m );
  double sum_lengths_p, sum_lengths_m;

  pt.uncovered_segments(piercings_p, depth_p, A, B, segments_p, sum_lengths_p, min_length);
  pt.uncovered_segments(piercings_m, depth_m, A, B, segments_m, sum_lengths_m, min_length);

  const double sum_lengths = (sum_lengths_p + sum_lengths_m);

  // if no uncovered segments of sufficient length
  if ( sum_lengths == 0. )
    return false;

  // pick interval uniform by length
  double t = Random::random_instance.generate_a_random_number();
  double fp = sum_lengths_p / sum_lengths;
  if (t < fp)
  {
    t /= fp;
    pt.pick_interval_uniformly( t, segments_p, sum_lengths_p, seg_start, seg_end );
  }
  else
  {
    t -= fp;
    t *= sum_lengths / sum_lengths_m;
    pt.pick_interval_uniformly( t, segments_m, sum_lengths_m, seg_start, seg_end );
    multiply(u, -1);
  }
  return true;
}

bool Spoke_Darts::generate_dart(double *dart, double *u, const double A, const double * c, Search_Structure *anchor_nbr,
                                const bool do_wheels, bool *dart_is_from_wheel_ptr, Wheel_Stats *wheel_stats,
                                bool *covered_sphere, const double radius_factor)
{

  assert( _rng->num_dim() == num_dim() );
  _rng->set_num_dim( num_dim() );

    // u = random dart direction
  _rng->sample_uniformly_from_unit_sphere_surface(u);

    // ==============================================
    // check if anchor point is covered by a sphere
    // optionally, wheel to find a new one that is uncovered
    // ==============================================
  bool dart_is_from_wheel = false;
  unsigned int pass = 0;
  bool valid_dart;

  do {
    // check if end of dart is inside the domain
    // transform direction u to a dart point
    // dart is point at distance A in direction u
    axpy(dart, A, u, c );
    bool in_domain = _domain.in_domain(dart);

    // achor covered?
    size_t covering_sphere;
    valid_dart = in_domain &&
      anchor_nbr->no_near_spheres( dart, covering_sphere, radius_factor );

    // anchor (or whole piercing) was covered,
    // try wheeling around the sphere
    bool covered_wheel = false;
    if (in_domain && !valid_dart && do_wheels)
    {
      if (wheel_stats)
        wheel_stats->start_clock();

      dart_is_from_wheel = true;

        // choose a new dart, and update u, from the wheel
        // assignment inside "if" is intended
      if ( (valid_dart = wheel( c, A, u, dart, anchor_nbr->sphere_container(), radius_factor )) )
        valid_dart = _domain.in_domain(dart);

        // in 2-d, if the spin was covered, then we're done
      if (!valid_dart)
      {
        covered_wheel = true;
        if (_num_dim == 2)
          *covered_sphere = true;
      }

      if (wheel_stats)
        wheel_stats->collect_stats(valid_dart);
    }

            // debug
//      plot_vertices_2d(mylibs._num_inserted_spheres, mylibs._spheres, 0, 0, 0, 0, xmin, xmax, iactive, false, 0, dart);

            // If dart was not valid, try the opposite diretion.
            //   This second-chance increases the chance of a hit above the proof criteria.
            //   However, it doesn't double it, so we can't count it as a miss.
    if (!valid_dart)
    {
      if (covered_wheel)
        pass=2;
      else
        multiply(u, -1.);
    }

    pass++;
  } while (!valid_dart && pass < 2);
    // end find uncovered anchor point

  if (dart_is_from_wheel_ptr) 
    *dart_is_from_wheel_ptr = dart_is_from_wheel;

  return valid_dart;
}


double Spoke_Darts::sample_from_spoke(Spoke_Length &sl, double r, double sample_start_min, double sample_start_max, double mid_1, double mid_2, double top, double mid_frac )
{
  // top appears to be unused

  // exactly zero length indicates some problem in some algorithms, but not all
  if ( sl.A_1 == sl.A_2 )
    return sl.A_1;
  
  // random bottom of spoke
  double bot = sample_start_min * r;
  if ( sample_start_min < sample_start_max )
  {
    double t = Random::random_instance.generate_a_random_number();
    bot += t * (sample_start_max - sample_start_min) * r;
  }

  // default placement at anchor
  double dx = sl.A;
  
  // stay away from any trimmed ends
  // smooth transition to zero placement probabilities
  
  double botp = bot;
  const double min_dist = bot - sl.A_1;
  
  double m1p, m2p, topp, zt;
  
  // was the spoke trimmed by any disks?
  if (sl.trimmed_top())
  {
    // ramp distribution
    //     *
    //    * *
    //   *   *
    //  *     *
    // A1  m  A2, m is mid_frac point
    // z b  t zt
    
    topp = sl.A_2 - min_dist * r;
    zt = sl.A_2;
    assert( botp < topp ); // enforced by ignore short spokes
    m1p = (mid_frac * topp + (1. - mid_frac) * botp);
    m2p = m1p;
  }
  else // !trimmed top
  {
    // ramp - box - ramp distribution
    //     ****
    //    *    *
    //   *      *
    //  *        *
    // A1  m1 m2  top
    // z b       t zt
    //
    
    m1p = mid_1 * r;
    m2p = mid_2 * r;
    topp = sl.A_2;
    zt = topp;
  }
  const double zb = sl.A_1;
  assert( sl.A_1 <= botp);
  assert( botp <= m1p);
  assert( sl.A_1 <  m1p);
  assert( m1p <= m2p);
  assert( m2p <  topp);
  assert( topp <= sl.A_2);
  
  // speedup: precompute these for untrimmed spokes
  
  const double peak = pow( m1p - zb, num_dim() - 1. );
  
  // by swept volume
  
  // bottom
  const double m1pv = pow( m1p - zb, num_dim() ) / num_dim();
  const double botpv = pow( botp - zb, num_dim() ) / num_dim();
  const double p_bot = m1pv - botpv;
  
  // box
  const double p_box = peak * (m2p - m1p);
  
  // top
  const double top_scale =  peak / pow( zt - m2p , num_dim() - 1);
  const double m2pv  = top_scale * pow( zt - m2p , num_dim() ) / num_dim();
  const double toppv = top_scale * pow( zt - topp, num_dim() ) / num_dim();
  const double p_top = m2pv - toppv;
  
  // total
  const double p_total = (p_box + p_bot + p_top);
  // random point of this distribution
  double t = Random::random_instance.generate_a_random_number();
  double p = t * p_total;
  
  // select from ramp-box-ramp
  // bot
  if (p < p_bot)
  {
    //              dx = zero + pow( (p * num_dim() + A_1_d), 1. / num_dim() );
    dx = zb + pow( num_dim() * ( p + botpv ), 1. / num_dim() );
    assert( dx >= sl.A_1 ); // roundoff?
    assert( dx >= botp ); // roundoff?
    assert( dx <= m1p ); // roundoff?
  }
  else
  {
    p -= p_bot;
    // top
    if ( p < p_top )
    {
      //                dx = topp - pow( p * num_dim() * pow(topp - m2p, num_dim() - 1) / peak, 1. / num_dim());
      //                dx = zerotop - pow( num_dim() * (p + top_cut_volume), 1. / num_dim() );
      dx = zt - pow( num_dim() * ( p + toppv ) / top_scale, 1. / num_dim() );
      assert( dx >= m2p ); // roundoff?
      assert( dx <= topp ); // roundoff?
      assert( dx <= sl.A_2 ); // roundoff?
    }
    // box
    else
    {
      p -= p_top;
      assert( p <= p_box * (1 + 1.0e-6));
      assert( p >= -p_box * 1.0e-6);
      dx = m1p + (m2p - m1p) * p / p_box;
      assert( dx >= m1p - r * 1.e-6);
      assert( dx <= m2p + r * 1.e-6);
    }
  }
  
  assert( dx >= sl.A_1 );
  assert( dx >= botp );
  assert( dx <= topp );
  assert( dx <= sl.A_2 );

  return dx;
}