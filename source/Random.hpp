// Random.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef RANDOM_HPP
#define RANDOM_HPP

// see http://www.thecodingforums.com/threads/re-rngs-a-double-kiss.720512/


#include "Point_Tool.hpp"
class Spoke_Length;

class Random : public Point_Tool
{

public:

  // sets the seed
  void initiate_random_generator(unsigned long x = 1391722129);

  // u in [0,1]
  double generate_a_random_number();

  // calculate some statistics
  void get_mean_and_variance(std::vector<double> & data, double &mean, double &var);


  // picking from intervals
  // return some random value favoring the middle of the interval from A to B,
  double random_middle(double A, double B);
  // pick uniformly from A to B
  double random_uniform(double A, double B);


  // these next rely on _num_dim being set

  // pick uniformly by volume, for a sphere or sector from radius r1 to r2
  // here t may be passed in
  double random_by_volume( double r );
  double random_by_volume( double r_inner, double r_outer );
  // From two annuli of differing extent
  //  e.g. for a bowtie, a double ended spoke; for a single segment trimmed at both ends
  double random_by_volume( bool &from_1, double r1_inner, double r1_outer, double r2_inner, double r2_outer );

// todo: decide where to put this:
  //double Balloon_Darts::sample_from_spoke_favored2(double r, Spoke_Length &sl, double bot, double mid_1, double mid_2, double mid_frac, double top)

  // s is a point, modified
  void sample_uniformly_from_box( double *s, const double *xmax, const double *xmin );

  // random point inside a ball
  void sample_uniformly_from_unit_sphere(double* dart);
  // random unit vector
  void sample_uniformly_from_unit_sphere_surface(double* dart);

  // find a random unit vector v that is orthogonal to a prior unit vector u
  void random_orthonormal_vector(const double *u, double *v);

  // for a spoke sl that goes through zero, a butterfly, pick one of the wings uniformly by area
  // returns true if we used the minus side, false if we used the plus side of sl
  bool random_spoke_side(const Spoke_Length *sl, Spoke_Length *onesided_sl);

  Random() : Point_Tool( Point_Tool::bad_num_dim() ) { initiate_random_generator(); }
  Random(size_t num_dim) : Point_Tool(num_dim) { initiate_random_generator(); }

// provided one a global, so we can set the seed just once
// but this is not required
  static Random random_instance;

private:
  // variables for Random number generator
	double Q[1220];
	int indx;
	double cc;
	double c; /* current CSWB */
	double zc;	/* current SWB `borrow` */
	double zx;	/* SWB seed1 */
	double zy;	/* SWB seed2 */
//	size_t qlen;/* length of Q array */
  
};



#endif