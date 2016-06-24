// Quality_Tool.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef QUALITY_TOOL_HPP
#define QUALITY_TOOL_HPP

#include <vector>
#include <math.h>
#include <limits.h>
#include <assert.h>

class Spheres;
class Domain;
class Search_Structure;


// public interface is just the data
// the public interface to the method to build it is in Quality_Tool
class Histogram
{
public:
// protected:
  std::vector<double> bins;
  unsigned int num_windowed_points, num_distances, num_histogram_distances, bins_per_rx;
  size_t dim, num_spheres;
  double window_frame, rx, bin_average, maxbin;
  double min_distance;
  bool look_periodic;
  
public:
  // constructor
  Histogram() :
  num_windowed_points(0), num_distances(0), num_histogram_distances(0), bins_per_rx(20), dim(2), num_spheres(0),
  window_frame(0.), rx(0.), bin_average(0.), maxbin(0.),
  look_periodic(false)
  {}

protected:
  friend class Quality_Tool;
  
  // return true if an error, such as too close of points
  bool build(double radius_factor, Spheres *spheres, Domain *domain, Search_Structure *search );
};



class Quality_Tool
{

public:    

  Quality_Tool() {}
  
  // radius_factor = 3 means report radii between 1 and 4=(1+3).
  bool build_histogram(Histogram &histogram, Spheres *spheres, Domain *domain, Search_Structure *search,  double radius_factor = 3.  )
  { return histogram.build( radius_factor, spheres, domain, search ); }
  
};

#endif