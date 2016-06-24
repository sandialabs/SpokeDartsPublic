// Random.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Random.hpp"
#include "Point_Tool.hpp"
#include "Spoke_Length.hpp"

Random Random::random_instance;



/////////////////////////////////////////////////////////////////////
// Random Methods
/////////////////////////////////////////////////////////////////////

void Random::initiate_random_generator(unsigned long x)
{
  #pragma region initiate The Random generator:
  //assert(sizeof (double) >= 54) ;

  cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
  size_t i;
  size_t qlen = indx = sizeof Q / sizeof Q[0];
  for (i = 0; i < qlen; i++) Q[i] = 0;

//  double c = 0.0, zc = 0.0, /* current CSWB and SWB `borrow` */ //zzyk - is this a bug? should we be using the global c and zc and zx and zy?
//  zx = 5212886298506819.0 / 9007199254740992.0, /* SWB seed1 */
//  zy = 2020898595989513.0 / 9007199254740992.0; /* SWB seed2 */
  c = 0.0, zc = 0.0, /* current CSWB and SWB `borrow` */
  zx = 5212886298506819.0 / 9007199254740992.0, /* SWB seed1 */
  zy = 2020898595989513.0 / 9007199254740992.0; /* SWB seed2 */
  
  size_t j;
  double s, t;   /* Choose 32 bits for x, 32 for y */
  if (x == 0) x = 123456789; /* default seeds */
  unsigned long y = 362436069; /* default seeds */

  /* Next, seed each Q[i], one bit at a time, */
  for (i = 0; i < qlen; i++)
  { /* using 9th bit from Cong+Xorshift */
    s = 0.0;
    t = 1.0;
    for (j = 0; j < 52; j++)
    {
      t = 0.5 * t; /* make t=.5/2^j */
      x = 69069 * x + 123;
//      y ^= (y << 13);
//      y ^= (y >> 17);
//      y ^= (y << 5);
      y = (y ^ (y << 13)) & 0xFFFFFFFFUL;
      y = (y ^ (y >> 17)) & 0xFFFFFFFFUL;
      y = (y ^ (y << 5)) & 0xFFFFFFFFUL;
      if (((x + y) >> 23) & 1) s = s + t; /* change bit of s, maybe */
    }  /* end j loop */
    Q[i] = s;
  } /* end i seed loop, Now generate 10^9 dUNI's: */
  #pragma endregion
}

double Random::generate_a_random_number()
{ 
  #pragma region generate a random number:
  /* Takes 14 nanosecs, Intel Q6600,2.40GHz */
  int i, j;
  double t; /* t: first temp, then next CSWB value */
  /* First get zy as next lag-2 SWB */
  t = zx - zy - zc;
  zx = zy;
  if (t < 0)
  {
    zy = t + 1.0;
    zc = cc;
  }
  else
  {
    zy = t;
    zc = 0.0;
  }

  /* Then get t as the next lag-1220 CSWB value */
  if (indx < 1220)
    t = Q[indx++];
  else
  { /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
    for (i = 0; i < 1220; i++)
    {
      j = (i < 30) ? i + 1190 : i - 30;
      t = Q[j] - Q[i] + c; /* Get next CSWB element */
      if (t > 0)
      {
        t = t - cc;
        c = cc;
      }
      else
      {
        t = t - cc + 1.0;
        c = 0.0;
      }
      Q[i] = t;
    }  /* end i loop */
    indx = 1;
    t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
  } /* end else segment; return t-zy mod 1 */
  
  return ((t < zy) ? 1.0 + (t - zy) : t - zy);  
  #pragma endregion
}

void Random::sample_uniformly_from_box( double *s, const double *xmax, const double *xmin )
{
  for (size_t d = 0; d < num_dim(); ++d)
    s[d] = xmin[d] + (xmax[d] - xmin[d]) * generate_a_random_number();
}


void Random::sample_uniformly_from_unit_sphere_surface(double* dart)
{
  assert(valid_num_dim());
  // unbiased method
  double sf = 0.0;
  for (size_t idim = 0; idim < num_dim(); idim++)
  {
    double sum(0.0);
    // select 12 random numbers from 0.0 to 1.0
    for (size_t i = 0; i < 12; i++) sum += generate_a_random_number();
    sum -= 6.0;
    dart[idim] = sum;
    sf += dart[idim] * dart[idim];
  }
  sf = 1.0 / sqrt(sf);
  for (size_t idim = 0; idim < num_dim(); idim++) 
    dart[idim] *= sf;
}

void Random::sample_uniformly_from_unit_sphere(double* dart)
{
  sample_uniformly_from_unit_sphere_surface(dart);  
  // point is on sphere surface, now project to interior

  double u = pow(generate_a_random_number(), 1.0 / num_dim());
  for (size_t idim = 0; idim < num_dim(); idim++) 
    dart[idim] *= u;
}

void Random::get_mean_and_variance(std::vector<double> & data, double &mean, double &var)
{
  #pragma region Get mean and average of a given distribution:  
  size_t num_experiments(data.size());
  mean = 0.0;
  for (size_t iexp = 0; iexp < num_experiments; iexp++) mean+= data[iexp];
  mean /= num_experiments;

  var = 0.0;
  for (size_t iexp = 0; iexp < num_experiments; iexp++)
  {
    double dvar = data[iexp] - mean;
    var += dvar * dvar;
  }
  var = sqrt(var/num_experiments);
  #pragma endregion
};


// find a random unit vector v that is orthogonal to unit vector u
void Random::random_orthonormal_vector(const double *u, double *v)
{
  assert( num_dim() < size_t(-1) );
  
  double dot(0.);
  double *w = new_point();
  do {
    sample_uniformly_from_unit_sphere_surface(w);
    dot = dot_product(u,w);
  }
  while (fabs(dot) > 0.9);
  // make v perpendicular to u
  // v = w - dot u
  axpy(v, -dot, u, w);
  const double norm_v = sqrt(1. - dot * dot);
  multiply(v, 1. / norm_v);
  assert( fabs(norm_squared(v) - 1.) < 0.0001 );
  delete_point(w);
}  

 
double Random::random_middle( double A, double B)
{  
  const double delta = ( generate_a_random_number() + generate_a_random_number() + generate_a_random_number() ) * (B-A) / 3.;
  return A + delta;
}  

double Random::random_uniform(double A, double B)
{
  const double t = generate_a_random_number();
  double m = t * (B-A) + A;
  return m;
}

double Random::random_by_volume( double r )
{
  assert(valid_num_dim());
  const double t = generate_a_random_number();
  return inverse_relative_volume( t, r );
}

double Random::random_by_volume( double r_inner, double r_outer )
{
  assert(valid_num_dim());
  const double t = generate_a_random_number();
  return inverse_relative_volume( t, r_inner, r_outer );
}

double Random::random_by_volume(bool &from_1, double r1_inner, double r1_outer, double r2_inner, double r2_outer )
{
  assert(valid_num_dim());
  const double v1_in  = relative_volume( r1_inner );
  const double v1_out = relative_volume( r1_outer );
  const double v_1 = v1_out - v1_in;
  const double v2_in  = relative_volume( r1_inner );
  const double v2_out = relative_volume( r2_outer );
  const double v_2 = v2_out - v2_in;
  const double t = generate_a_random_number();
  const double p = t * (v_1 + v_2);
  double v = p + v1_in;
  from_1 = true;
  if ( p < v_1 )
  {
    v = p - v_1 + v2_in;
    from_1 = false;
  }
  const double r = pow( v, 1. / num_dim());
  return r;
}


   /*
     // d+1 triangle-distribution
    double g = generate_a_random_number();
    for (unsigned int dd = 0; dd < num_dim(); ++dd)
    {
      g += generate_a_random_number();
    }
    const double d_theta = g * best_interval_length / (num_dim()+1);
     */

bool Random::random_spoke_side(const Spoke_Length *sl, Spoke_Length *onesided_sl)
{
  assert(sl->A_1 <= 0.);
  assert(sl->A_2 >= 0.);
  const double len_1 = -sl->A_1;
  const double len_2 =  sl->A_2;
  const double area_1 = relative_volume(len_1); 
  const double area_2 = relative_volume(len_2); 
  const double t = generate_a_random_number();
  bool pick_minus = t  <  area_1 / (area_1 + area_2);
  
  // set onesided to the corresponding side
  if (pick_minus)
  {
    // replace all A_2 like quantities with -A_1
    // replace all A_1 quantities with zero
    onesided_sl->A = (sl->A < 0) ? fabs(sl->A) : 0.;
    onesided_sl->A_2 = fabs(sl->A_1);
    onesided_sl->A_2_init = fabs(sl->A_1_init);
    onesided_sl->A_2_max = fabs(sl->A_1_min);
    onesided_sl->A_1 = 0;
    onesided_sl->A_1_init = 0;
    onesided_sl->A_1_min = 0;
  }
  else
  {
    // replace A_1 with zero
    onesided_sl->A = (sl->A > 0) ? sl->A : 0.;
    onesided_sl->A_2 = sl->A_2;
    onesided_sl->A_2_init = sl->A_2_init;
    onesided_sl->A_2_max = sl->A_2_max;
    onesided_sl->A_1 = 0;
    onesided_sl->A_1_init = 0;
    onesided_sl->A_1_min = 0;
  }
  onesided_sl->_r = sl->_r;

  
  return pick_minus;
}

