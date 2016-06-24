// Favored_Spoke_MPS.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

class Favored_Spoke_MPS
{
public:
  // =========
  // INTERFACE
  // call this function to create an MPS distribution using spokes
  // =========

private:
  // minimum length of a spoke that we insert a point along
  double min_length();
};

double Spoke_MPS::min_length()
{
  double min_length = 0.;
  if (greedy_pass)
    min_length = B-A;
  else if (!soft_skip_short_spokes)
  {
    min_length = skip_short_spokes * r;
  }
  else if (soft_skip_short_spokes)
  {
    double st = generate_a_random_number();
    // zzyk. research question
    // Here st is cubed, regardless of dimension. It was chosen experimentally for d=2.
    // Perhaps the power should be dimension dependent, num_dim+1?
    min_length = (st * st * st * (soft_skip_max - soft_skip_min) + soft_skip_min) * r;
  }
  return min_length;
}

to pick a piercing segment

  Spoke_Darts sd( _spheres, _domain );
  
  const double *p = _spheres[ci];
  KD_Tree trim_spheres_tree; // extended neighbors
  size_t *trim_array = trim_spheres_tree.array();
  size_t &trim_size = trim_spheres_tree.size();

  pick_segment_piercing(const size_t ci, const double r, const double A, const double B, double * u,
                          double & seg_start, double & seg_end, 
                          trim_array, trim_size,
                          greedy_pass, min_length());



double Balloon_Darts::sample_from_spoke_favored2(double r, Spoke_Length &sl, double bot, double mid_1, double mid_2, double mid_frac, double top)
{
  double t = mylibs.generate_a_random_number(); // random parameter for placement
  double dx = sl.A; // default placement at anchor
  
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
  
  const double peak = pow( m1p - zb, _num_dim - 1. );
  
  // by swept volume
  
  // bottom
  const double m1pv = pow( m1p - zb, _num_dim ) / _num_dim;
  const double botpv = pow( botp - zb, _num_dim ) / _num_dim;
  const double p_bot = m1pv - botpv;
  
  // box
  const double p_box = peak * (m2p - m1p);
  
  // top
  const double top_scale =  peak / pow( zt - m2p , _num_dim - 1);
  const double m2pv  = top_scale * pow( zt - m2p , _num_dim ) / _num_dim;
  const double toppv = top_scale * pow( zt - topp, _num_dim ) / _num_dim;
  const double p_top = m2pv - toppv;
  
  // total
  const double p_total = (p_box + p_bot + p_top);
  double p = t * p_total;
  
  // select from ramp-box-ramp
  // bot
  if (p < p_bot)
  {
    //              dx = zero + pow( (p * _num_dim + A_1_d), 1. / _num_dim );
    dx = zb + pow( _num_dim * ( p + botpv ), 1. / _num_dim );
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
      //                dx = topp - pow( p * _num_dim * pow(topp - m2p, _num_dim - 1) / peak, 1. / _num_dim);
      //                dx = zerotop - pow( _num_dim * (p + top_cut_volume), 1. / _num_dim );
      dx = zt - pow( _num_dim * ( p + toppv ) / top_scale, 1. / _num_dim );
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
