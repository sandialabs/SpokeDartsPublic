// Spoke_Length.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// define the extent of spokes and trimmed spokes

#ifndef SPOKE_LENGTH
#define SPOKE_LENGTH

class Spoke_Length
{
public:

  // The spoke goes from distance A_1 to distance A_2, with anchor point A between them.
  // min and max are the longest initial spoke we might start with
  // init is the initial spoke length we picked, before trimming,
  // A_1 and A_2 are the distances after trimming

  // _min and _max are abstract parameters like 1 and 2
  // all other quantities are multiples of r
  double A, A_1_min, A_2_max, A_1_init, A_2_init, A_1, A_2, _r;
  
  Spoke_Length(double A_1_min_in = 1., double A_2_max_in = 2., double r = 1.)
  : A_1_min(A_1_min_in), A_2_max(A_2_max_in),
    A_1_init(A_1_min*r), A_2_init(A_2_max*r), A_1(A_1_init), A_2(A_2_init)
  {
    reset(r); // redundant safety
  }

  void reset(double r)
  {
    _r = r;
    A_1 = A_1_init = A_1_min*r;
    A_2 = A_2_init = A_2_max*r;
    A = r;
  }

  void reset(double A_1_min_in, double A_2_max_in, double r)
  {
    A_1_min = A_1_min_in;
    A_2_max = A_2_max_in;
    reset(r);
  }
  
  const double &r() const {return _r;}
  
  bool trimmed_bottom() const
  {
    return A_2 < A_2_init;
  }

  bool trimmed_top() const
  {
    return A_2 < A_2_init;
  }

  bool was_trimmed() const
  {
    return trimmed_top() || trimmed_bottom();
  }
  
};

#endif