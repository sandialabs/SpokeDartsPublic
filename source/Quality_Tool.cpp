// Quality_Tool.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Quality_Tool.hpp"


#include "Spheres.hpp"
#include "Domain.hpp"
#include "Search_Structure.hpp"
#include "Global_Container.hpp"
#include "Sphere_Container_Array.hpp"


bool Histogram::build(double radius_factor, Spheres *spheres, Domain *domain, Search_Structure *search )
{
  bool error(false);
  Point_Tool pt( domain->num_dim() );
  look_periodic = domain->is_periodic();
  
  dim = domain->num_dim();
  num_spheres = spheres->num_real_spheres();
  
  // only collect statistics for points inside the window_frame.
  // it should something more like 4 to eliminate all edge effects, but that often reduces the window to zero.
  rx = spheres->radius(0); // radius of first point, might not be constant

  bins_per_rx = 20;
  // max_h should be less than the window frame, or half the domain size if periodic
  const double max_h = radius_factor + 1./bins_per_rx;
  const double max_neighbor = (1. + max_h) * rx;
  const unsigned int numbins = (unsigned int)ceil( max_h * bins_per_rx );
  bins.assign(numbins,0.);
  
  // don't count distances from points outside the frame
  // the reason is to avoid boundary effects messing up the histogram
  // for non-periodic domains, don't trust histogram values larger than the frame
  // is_periodic ? 0 : 3 are the safest values, but there is often not enough data to use 3
  if (look_periodic)
    window_frame = 0.;
  else
  {
    window_frame = 3.;
    const double frame_max = 0.3 / rx;
    if ( window_frame > frame_max )
      window_frame = frame_max;
  }

  num_windowed_points = 0;
  num_distances = 0;
  num_histogram_distances = 0;
  
  // build histogram
  Global_Container outer(*spheres); // doesn't include ghosts
  Sphere_Container_Array neighbors(*spheres);
  
  min_distance = 1.;
  size_t num_errors=0;
  
  for (size_t i = outer.first(); i != outer.bad_sphere_index(); i = outer.next())
  {
    const double *p = (*spheres)[i];
    
    if ( look_periodic && !pt.in_box(p, domain->xmax(), domain->xmin(), rx * (1e-4) ) )
    {
      error = true;
      std::cerr << "Domain error, real disk i:" << i;
      pt.print_point( p, std::cerr );
      std::cerr << " is outside the domain box: ";
      pt.print_point(domain->xmin(), std::cerr);
      std::cerr << " x ";
      pt.print_point(domain->xmin(), std::cerr);
      std::cerr << std::endl;
    }
    
    if ( look_periodic || pt.in_box( p, domain->xmax(), domain->xmin(), window_frame ))
    {
      num_windowed_points++;
      // todo : speedup, rewrite all_near_spheres so it just returns the distances, rather than the ids.
      search->all_near_spheres( &neighbors, p, max_neighbor );
      
      for (size_t j = neighbors.first(); j != neighbors.bad_sphere_index(); j = neighbors.next())
      {
        const double *q = (*spheres)[j];
        
        double h = sqrt( pt.distance_squared( p, q ) );
        
        if ( h < min_distance )
          min_distance = h;
        
        if ( h < rx)
        {
//          if ( h < rx * (1. - 1e-3) - 1e-7 )
          if ( h < rx * (1. - 1e-2) - 1e-7 )
          {
            error = true;
            if ( ++num_errors < 10)
            {
              std::cerr << "Conflict distance error, disks i:" << i << " j:" << j << " at distance:" << h << " < rx:" << rx << ", distance = " << h / rx << " r." << std::endl;
              std::cerr << "  " << i << " "; pt.print_point( p, std::cerr ); std::cerr << std::endl;
              std::cerr << "  " << j << " "; pt.print_point( q, std::cerr ); std::cerr << std::endl;
            }
          }
          h = rx; // so it gets counted as a close distance
        }
        num_distances++;
        unsigned int k = (unsigned int) ( floor( ( h/rx - 1. ) * bins_per_rx ) ); // floor
        if (k < numbins)
        {
          bins[k]++;
          num_histogram_distances++;
        }
      }
    }
  }
  min_distance = min_distance / rx;
  
  // normalize histogram
  maxbin=0.;
  for (unsigned int k = 0; k < numbins; k++)
  {
    // by surface area
    // std::cout << " scale:" << pow( (k + bins_per_rx + 0.5), dim - 1.) << "  ";
    bins[k] /= pow( (k + bins_per_rx + 0.5), pt.num_dim() - 1.);
    
    // find max value for lineaer scaling of display
    if (bins[k] > maxbin)
      maxbin = bins[k];
  }
  
  // bin average
  bin_average = 0.;
  unsigned int num_k = 0;
  for (unsigned int k = 0; k < numbins; k++)
  {
    if ( k > bins_per_rx )
    {
      num_k++;
      bin_average += bins[k];
    }
  }
  bin_average /= num_k;
  
  return error;
}


