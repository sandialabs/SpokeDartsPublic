// Ghost_Spheres.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Ghost_Spheres.hpp"
#include "Sphere_Container_Array.hpp"

Ghost_Spheres::Ghost_Spheres(size_t max_num_spheres, Domain &domain)
: Spheres(max_num_spheres, domain), _frame(0.),
_two_power_dim ( ( size_t(1))  << domain.num_dim() ), _periodic_copies( domain.num_dim(), _two_power_dim),
_permanent_ghosts( true ), _num_ghosts(0),
_periodic_copies_size(0)
{
  _periodic_copies.new_memory( _two_power_dim );
  _added_spheres = new size_t[ _two_power_dim ];
}
Ghost_Spheres::Ghost_Spheres(Domain &domain)
: Spheres(domain), _frame(0.),
_two_power_dim ( ( size_t(1))  << domain.num_dim() ), _periodic_copies(domain.num_dim(), _two_power_dim),
_permanent_ghosts( true ), _num_ghosts(0),
_periodic_copies_size(0)
{
  _periodic_copies.new_memory( _two_power_dim );
  _added_spheres = new size_t[ _two_power_dim ];
}

Ghost_Spheres::~Ghost_Spheres()
{
  // _periodic_copies will be freed by its own destructors
  dematerialize_ghosts();
  delete [] _added_spheres;
}

void Ghost_Spheres::dematerialize_ghosts()
{
  for ( size_t i = first_ghost_index(); i < max_num_spheres(); ++i)
  {
    _pt.delete_sphere( _sphere_array[i] );
    // or throw back to the pool;
  }
  _num_ghosts = 0;
}

void Ghost_Spheres::make_periodic_point_copies( const double *c )
{
  _periodic_copies_size = 0;
  if ( no_ghosts_needed( c ) )
  {
    _pt.assign(_periodic_copies[0], c); // no memory allocation needed
    _periodic_copies_size = 1;
    _real_index = 0;
    _orig_index = 0;
    return;
  }
  double *c_copy = _pt.new_sphere();
  _pt.assign(c_copy, c);
  _pt.radius(c_copy) = 0.; assert(0); // currently we shouldn't be calling this. Note the radius is bad.
  make_periodic_copies_recurse( 0, true, true, c_copy );
}

void Ghost_Spheres::make_periodic_sphere_copies( const double *s )
{
  _periodic_copies_size = 0;
  if ( no_ghosts_needed( s ) )
  {
    _pt.assign_sphere(_periodic_copies[0], s); // no memory allocation needed
    _periodic_copies_size = 1;
    _real_index = 0;
    _orig_index = 0;
    return;
  }
  double *s_copy = _pt.new_copy_sphere(s);
  make_periodic_copies_recurse( 0, true, true, s_copy );
  _pt.delete_sphere(s_copy);
}


void Ghost_Spheres::make_periodic_copies_recurse( size_t d, bool is_orig, bool is_real, double *c )
{
  // fathom
  if (d == num_dim())
  {
    assert( _periodic_copies_size < _two_power_dim ); // array bounds check
    // for very large radii, we might need more than one copy in each direction! 3^d or more!
    _pt.assign_sphere(_periodic_copies[_periodic_copies_size], c); // no memory allocation needed
    if (is_real)
      _real_index = _periodic_copies_size;
    if (is_orig)
      _orig_index = _periodic_copies_size;
    ++_periodic_copies_size;
  }
  // recurse
  else
  {
    const double period = domain().xmax()[d] - domain().xmin()[d];

    // guarantee exactly one permutation is considered real
    size_t real_index = 2;
    if (is_real)
    {
      if ( c[d] - period >= domain().xmin()[d] )
        real_index = 0;
      else if ( c[d] >= domain().xmin()[d] )
        real_index = 1;
    }

    // orig
    if (in_frame(d, c[d]))
      make_periodic_copies_recurse(
        d+1,
        is_orig,
        is_real && real_index == 1,
        c);

    // orig - period
    const double saved_coordinate = c[d];
    c[d] -= period;
    if (in_frame(d, c[d]))
      make_periodic_copies_recurse(
        d+1,
        false,
        is_real && real_index == 0, 
        c);

    // orig + period
    c[d] = saved_coordinate + period;
    if (in_frame(d, c[d]))
      make_periodic_copies_recurse(
        d+1,
        false,
        is_real && real_index == 2, 
        c);
    c[d] = saved_coordinate;
  }
}

size_t Ghost_Spheres::add_permanent_ghosts(const double *c, size_t &real_sphere_index )
{
  size_t orig_sphere_index( bad_sphere_index() );
  make_periodic_sphere_copies( c );
  if (_periodic_copies_size)
  {
    for ( size_t i = 0; i < _periodic_copies_size; ++i )
    {
      size_t spheres_index;
      if (i == _real_index)
        spheres_index = real_sphere_index = add_real_sphere_to_array( _periodic_copies[i] );
      else
        spheres_index = add_ghost_sphere_to_array( _periodic_copies[i] );
      if (i == _orig_index)
        orig_sphere_index = spheres_index;
      // add them all
      //    if ( (i != _real_index) && (i != _orig_index) )
      _added_spheres[ _num_added_spheres++ ] = spheres_index;
    }
  }
  else
  {
    orig_sphere_index = real_sphere_index = add_real_sphere_to_array( c );
  }
  return orig_sphere_index;
}


bool Ghost_Spheres::ghost_to_real(const double *g, double *c )
{
  // this used to check that only one addition or subtraction is needed
  // but this isn't required anymore!
  bool changed = false;
  for (size_t d = 0; d<num_dim(); ++d)
  {
    const double period = domain().xmax()[d] - domain().xmin()[d];
    c[d] = g[d];
    while (c[d] > domain().xmax()[d])
    {
      c[d] -= period;
      changed = true;
    }
    while (c[d] < domain().xmin()[d])
    {
      c[d] += period;
      changed = true;
    }
    assert( c[d] >= domain().xmin()[d] );
    assert( c[d] <= domain().xmin()[d] + period + 1e-15); // <= xmax
  }
  return changed;
}

size_t Ghost_Spheres::add_sphere(const double *sphere, size_t &real_index)
{
  _num_added_spheres = 0;
  if (!domain().is_periodic())
  {
    return Spheres::add_sphere(sphere, real_index);
  }
  if (_permanent_ghosts)
  {
    return add_permanent_ghosts( sphere, real_index );
  }
  // temporary_ghosts, add the real copy and the original only
  {
    double *s = _pt.new_sphere();
    _pt.radius(s) = _pt.radius(sphere);
    bool changed = ghost_to_real(sphere, s);
    real_index = Spheres::add_real_sphere_to_array(s);
    // create a temporary ghost in the array, if needed
    return ( !changed ? real_index : add_ghost_sphere_to_array( sphere ) );
  }
}
