// Grid_Container.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// This is the data for the grid, basic add item 

#ifndef GRID_CONTAINER_HPP
#define GRID_CONTAINER_HPP

#include "Sphere_Container.hpp"
#include "Point_Tool.hpp"
#include "Sphere_Container_Array.hpp"

class Grid_Container : public Sphere_Container
{

public:

  // Constructor
  // search_distance is a lower bound on the cell side length
  // todo, extend to a version where each dimension has a different range
  Grid_Container(Spheres &spheres, double search_distance, double xmax, double xmin);

  virtual ~Grid_Container()
  { delete_memory(); }
  
  // Contains the indices of the spheres in the cell
  typedef std::vector< size_t > Cell;
  // for retrieving a particular cell
  // array of indices, of length num_dim
  typedef int *Cell_Index;

  // Geometric cell goes from index*_cell_side to (index+1)*_cell_side
  
  // return the cell (that would be) containing point _spheres[s]
  Cell *cell( size_t isphere );

  // return the cell (that would be) containing point p
  Cell *cell( double *p );

  // Sphere_Container interface functions
  
  // get rid of contents
  virtual void clear();

  // add one sphere
  virtual void add_sphere(size_t isphere);
  
  // number of elements
  virtual size_t size() const {return _array.size();}

  // access the array form of the grid
  const Sphere_Container_Array &array() { return _array; }

  // grid_container[i]
  //
  // rvalue access only through operator []
  // size_t sphere_index = grid_container[9];
  const size_t & operator[](size_t i) const
  {
    assert(i < size());
    return _array[i];
  };
  // lvalue 
  size_t & operator[](std::size_t i)
  {
    return _array[i];      /* actual access, e.g. return mVector[idx]; */
  };
  

protected:
  // methods

  // add node to grid
  // the node has to be added to the array separately
  void add_sphere_to_grid(const size_t isphere);
  
  // memory allocation
  void new_memory();
  void delete_memory();

protected:

  // data

  // spacing parameters passed in
  double _cell_side; // rounded from input
  double _xmax, _xmin;
  
  // derived spacing parameters
  double _cell_diagonal;

  // number of cells in each dimension,
  // i.e. _cell_side * _num_division = _xmax-x_min
  size_t _axis_divisions;
  size_t _num_cells; // _axis_divisions^num_dim

  // grid, multi-dim array of vectors
  Cell *_cells;
  
  // handy cell index
  // reused to avoid frequent memory allocations
  Cell_Index _ci;
  
  // lookup a cell
  // idiom, for point p
  // Cell_Index ci = new int[ num_dim() ];
  // point_to_cell_index(p, ci);
  // Cell * = & _cells[ cell_index_to_linear_index(ci) );
  // delete [] ci;
  int coordinate_to_axis_index( double p ) const;
  void point_to_cell_index( const double *p, Cell_Index ci) const; // ci is changed
  size_t cell_index_to_linear_index( const Cell_Index ci) const; // ci is const
  Cell *cell_from_cell_index( const Cell_Index ci ) const;
  Cell *cell_from_linear_index( size_t li ) const;
  
  double upper_coordinate( Cell_Index ci, size_t d )
  {
    return (ci[d] + 1) * _cell_side + _xmin;
  }
  double lower_coordinate( Cell_Index ci, size_t d )
  {
    return ci[d] * _cell_side + _xmin;
  }

  
  // Search_Grid
  // iterate through the neighbors
  // void next_neighbor( Cell_Index center, size_t &state, Cell_Index next );
  
  // straight array of all the indices in the container
  // takes space and time to build, but fast to access
  Sphere_Container_Array _array;
  
  virtual void reset_iterator_it( size_t &it ) const
  {
    return _array.reset_iterator_it( it );
  }
  virtual size_t next_it( size_t &it ) const
  {
    return _array.next_it( it );
  }

private:

  void clear_grid();

};

inline
int Grid_Container::coordinate_to_axis_index( double p ) const
{
  int i = (int) ( floor( (p - _xmin) / _cell_side )  );
  // if we want to enforce that p is in the box, then
  //  if ( i == _axis_divisions )
  //    i = _axis_divisions - 1;
  // assert( i < _axis_divisions ); // fails if p < _xmin or p >= xmax + _cell_side
  return i;
}

inline
Grid_Container::Cell * Grid_Container::cell( size_t isphere )
{
  assert( _spheres.is_valid_sphere(isphere) );
  const double *p = _spheres[isphere];
  assert(p);
  point_to_cell_index( p, _ci);
  size_t li = cell_index_to_linear_index(_ci);
  Cell * c = & _cells[ li ];
  return c;
}

inline
Grid_Container::Cell * Grid_Container::cell( double *p )
{
  point_to_cell_index( p, _ci);
  size_t li = cell_index_to_linear_index(_ci);
  assert( li >=0 && li < _num_cells );
  Cell * c = & _cells[ li ];
  return c;
}

inline
Grid_Container::Cell *Grid_Container::cell_from_cell_index( const Cell_Index ci ) const
{
  size_t li = cell_index_to_linear_index(ci);
  Cell * c = & _cells[ li ];
  return c;
}

inline
Grid_Container::Cell *Grid_Container::cell_from_linear_index( size_t li ) const
{
  Cell * c = & _cells[ li ];
  return c;
}



#endif
