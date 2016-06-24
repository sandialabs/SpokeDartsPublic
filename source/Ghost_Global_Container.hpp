// Ghost_Global_Container.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef GHOST_GLOBAL_CONTAINER
#define GHOST_GLOBAL_CONTAINER

#include "Global_Container.hpp"
#include "Ghost_Spheres.hpp"

class Ghost_Global_Container: public Global_Container
{
public:
  // create a ghost or non-ghost global container
  static Global_Container *new_global_container( Spheres& spheres );

  Ghost_Global_Container(Ghost_Spheres &ghost_spheres) : Global_Container(ghost_spheres), _ghost_spheres(ghost_spheres) {}

  // mostly the same as with non-ghosts, except the following
  
  virtual size_t size() const 
  { return _ghost_spheres.size(); }
  
  virtual void set( size_t i )
  {
    if ( i < _ghost_spheres.num_real_spheres() )
      _it = i;
    else if ( i < _ghost_spheres.first_ghost_index() )
      _it = _ghost_spheres.first_ghost_index();
    else
      _it = i;
  }


protected:

  // return the current value of it, and it is advanced to the next valid value or to bad_sphere_index
  virtual size_t next_it( size_t &it ) const
  {
    // go over all the real spheres, then over all the ghosts
    // as with most iterators, changing the contents of the _ghosts_spheres may invalidate the iterator

    // bad_sphere_index, unchanging
    if (it == _spheres.bad_sphere_index())
      return it;

    // real
    size_t ret_it = it++;
    // if ( it < _ghost_spheres.Spheres::size() ) , done
    
    // first ghost
    if ( it == _ghost_spheres.Spheres::size() )
    {
      it = _ghost_spheres.first_ghost_index();
    }
    // valid ghost
    if ( it == _ghost_spheres.max_num_spheres())
    {
      it = _spheres.bad_sphere_index();
    }

    return ret_it;
  }


private:
  // _spheres is also available
  Ghost_Spheres &_ghost_spheres;

};

inline
Global_Container *Ghost_Global_Container::new_global_container( Spheres& spheres )
{
  return (spheres.is_ghost() ? new Ghost_Global_Container( *(spheres.cast_to_ghost() ) ) : new Global_Container(spheres));
}


#endif
