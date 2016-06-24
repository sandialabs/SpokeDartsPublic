 // Ghost_Spheres.hpp

#ifndef GHOST_SPHERES_HPP
#define GHOST_SPHERES_HPP

// ghost spheres for periodic domains
// two ways: temporary ghosts, and permanent ones

#include <vector>
#include "Spheres.hpp"
#include "Domain.hpp"

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

class Frame;

class Ghost_Spheres : public Spheres
{
public:
  // ================================================
  // Setup
  // These settings affect the implementation,
  // should be set only once
  // ================================================

  Ghost_Spheres(size_t max_num_spheres, Domain &domain);
  Ghost_Spheres(Domain &domain);
  virtual ~Ghost_Spheres();

  void set_frame( double frame ) {_frame=frame;}
  double frame() const {return _frame;}
  void set_permanent_ghosts( bool permanent ) {_permanent_ghosts = permanent;}

  bool is_permanent_ghosts() const { return domain().is_periodic() && _permanent_ghosts; }
  bool is_temporary_ghosts() const { return domain().is_periodic() && !_permanent_ghosts; }

  // ================================================
  // overloaded from Spheres, takes care of ghosting
  // For a sphere that may or may not be in the (non-periodic) domain,
  // Add it, and any necessary ghost copies of it, to the Spheres
  // Return the index of the original point, and the real point, in Spheres; these may be the same.
  // Also "return" in the workspace all of the added spheres indices
  virtual size_t add_sphere(const double *sphere, size_t &real_index);
  size_t *_added_spheres;
  size_t _num_added_spheres;
  
  virtual bool remove_sphere(size_t &i)
  {
    if (i == bad_sphere_index())
      return false;

    // remove real sphere
    if (i < Spheres::size())
    {
      bool ret_val = Spheres::remove_sphere(i); // i is advanced
      if ( i == bad_sphere_index() && _num_ghosts > 0 )
        i = first_ghost_index();
      return ret_val;
    }
    
    if (i >= first_ghost_index() && i < _max_num_spheres )
    {
      // replace this sphere with the top-indexed ghost
      double *removed_sphere = _sphere_array[i];
      _sphere_array[i] = _sphere_array[first_ghost_index()];
      _pt.delete_sphere(removed_sphere);
      // forget the top-indexed ghost
      if (_num_ghosts>0)
        --_num_ghosts;
      // increment i, because it was replaced by a prior sphere
      if (++i == max_num_spheres())
        i = bad_sphere_index();
      return true;
    }
    return false;
  }


  virtual size_t size() const {return Spheres::size() + _num_ghosts;}

  // class specific
  // index of the first ghost in the _sphere_array, or _max_num_spheres if none (not bad_index)
  size_t first_ghost_index() const {return _max_num_spheres - _num_ghosts;} 
  // index of the first ghost in the _sphere_array
  size_t num_ghosts() const {return _num_ghosts;} 

  // ================================================


  //=================================================
  // temporary ghosts
  //=================================================

  // Just like Search_Structure definition.
  // It materializes nearby ghosts if necessary, creating those within distance dist to c
  // These are temporary ghosts.
  // The search is used to find/create these points.
  // The nearby ghosts, and nearby real points, are returned within the sub_search Search_Structure.

  // get rid of all temporary ghosts in the Spheres datastructure
  void dematerialize_temporary_ghosts()
  { if (domain().is_periodic() && !_permanent_ghosts) dematerialize_ghosts(); }
  // get rid of all materialized ghosts, whether temporary or not. Use caution before calling, as this is just for cleanup.
  void dematerialize_ghosts();

  // make all the periodic copies of c that are within the frame
  // the results are "returned" by setting the publically accessible period_copies array of coordinates
  // and the size, and which index is the original and real copy
  void make_periodic_point_copies( const double *c );
  void make_periodic_sphere_copies( const double *s );
  const Sphere_Array &periodic_copies() const { return _periodic_copies; }
  const size_t &periodic_copies_size() const { return _periodic_copies_size; }
  const size_t orig_index() const {return _orig_index;}
  const size_t real_index() const {return _real_index;}
  const bool is_orig_index(size_t i) const {return _orig_index == i;}
  const bool is_real_index(size_t i) const {return _real_index == i;}

  // for a sphere g that might be outside the domain, convert it into a real point, using the periodicity of the domain
  // return true if a real point was created, otherwise left as is
  // The real one is the copy that has the smallest postion that is still bigger than the domain min
  bool ghost_to_real(const double *g, double *c);

  virtual bool is_valid_sphere( size_t i ) const
  {
    return (
            // valid real sphere index
            Spheres::is_valid_sphere(i) ||
            // valid ghost sphere index
            (i >= first_ghost_index() && i < max_num_spheres())
            );
  }
  
  
  // true if the d'th coordinate x is within the domain + frame
  // used by Search_Tree and the like
  bool in_frame( size_t d, double x ) const
  { return x >= domain().xmin()[d] - _frame && x < _domain.xmax()[d] + _frame; }
  // true if all coordinates of the point is in the frame
  bool in_frame( const double *x ) const
  {
    for (size_t d = 0; d < num_dim(); ++d)
      if (!in_frame(d, x[d]))
        return false;
    return true;
  }
  // true if the d'th coordinate x is within the domain - frame,
  // and hence any ghost of it will be outside the frame
  // used by Search_Tree and the like
  bool no_ghosts_needed( const double *x ) const
  { return no_ghosts_needed( x, _frame ); }
  bool no_ghosts_needed( const double *x, double frame_size ) const
  {
    for (size_t d = 0; d < num_dim(); ++d)
    {
      if ( x[d] < domain().xmin()[d] + frame_size || 
           x[d] > domain().xmax()[d] - frame_size )
        return false;
    }
    return true;
  }

  // For use by Ghost_Global_Container
  // is this a Spheres or a Ghost_Spheres?
  virtual bool is_ghost() const
  {return true;}
  // is this particular spheres a ghost?
  virtual bool is_ghost(size_t i) const
  {return i != bad_sphere_index() && i >= first_ghost_index();}
  virtual Ghost_Spheres* cast_to_ghost() const {return const_cast<Ghost_Spheres*>(this);}
  virtual const Ghost_Spheres* cast_to_const_ghost() const {return this;}

protected:
  friend class Search_Ghost;
  friend class Ghost_Global_Container;
  friend class Global_Container;
  friend class Search_Ghost_Array;

  // puts the sphere in the array of ghosts; the caller must ensure that it is really a ghost
  size_t add_ghost_sphere_to_array(const double *s )
  {
    _num_ghosts++;
    size_t new_index = first_ghost_index();
    assert( new_index > Spheres::size() ); // memory exhausted after this one

    assert( domain().out_bounds_tol( s ) ); // bad position, point is inside the domain box

    _sphere_array[ new_index ] = _pt.new_copy_sphere(s);
    return new_index;
  }
  
protected:

  // the size of the box around the domain where ghosts will be created.
  // farther ghosts will be ignored
  bool _permanent_ghosts;
  double _frame;
  // keep _two_power_dim before _periodic_copies
  size_t _two_power_dim; // two to the power n, the number of periodic copies a point could have
  size_t _periodic_copies_size;
  Sphere_Array _periodic_copies;   // allocated array of spheres, periodic copies of a point
  // index into _periodic_copies of the coordinates passed in, and the one in the domain
  size_t _orig_index, _real_index;

  // number of ghosts permanently stored at the end of the spheres array
  size_t _num_ghosts;

  // add the periodic copies of c to the Spheres
  // return the index into Spheres of the original coordinates c, and the real copy too
  size_t add_permanent_ghosts( const double *c, size_t &real_sphere_index );

private:
  // c really does get changed
  void make_periodic_copies_recurse( size_t d, bool is_orig, bool is_real, double *c );
};

#endif