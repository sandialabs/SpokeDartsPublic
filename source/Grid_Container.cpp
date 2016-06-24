// Grid_Container.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Grid_Container.hpp"
#include <assert.h>

Grid_Container::Grid_Container(Spheres &spheres, double search_distance, double xmax, double xmin)
: Sphere_Container( spheres ), _array( spheres ), _xmax(xmax), _xmin(xmin), _cells(0), _ci(0)
{
  assert( _xmax > _xmin );
  const double axis_length = _xmax - _xmin;
  assert( axis_length >= 0 );
  _axis_divisions = (size_t) ceil( axis_length / search_distance );
   // _axis_divisions^num_dim
  _num_cells = (size_t) ceil( pow( _axis_divisions, num_dim() ));
  _cell_side = axis_length / _axis_divisions;
  _cell_diagonal = sqrt(num_dim()) * _cell_side;

  new_memory();
  
};

void Grid_Container::new_memory()
{
  // delete old memory, if any
  delete_memory();

  // allocate new
  _cells = new Cell[_num_cells];
  _ci = new int[num_dim()];
}

void Grid_Container::delete_memory()
{
  delete []  _cells;  _cells = 0;
  delete [] _ci;      _ci = 0;
}

void Grid_Container::clear()

{
  clear_grid();

  _array.clear();
}
void Grid_Container::clear_grid()
{
  // traverse the grid and empty the vectors
  for (size_t li = 0; li < _num_cells; ++li )
    _cells[li].clear();
}

void Grid_Container::add_sphere(const size_t isphere)
{  
  // add to array
  _array.add_sphere(isphere);

  add_sphere_to_grid( isphere );
}

void Grid_Container::add_sphere_to_grid(const size_t isphere)
{
  // grid

  assert( _spheres.is_valid_sphere( isphere ) );

  Cell *c = cell(isphere);
  c->push_back( isphere );
}


void Grid_Container::point_to_cell_index( const double *p, Cell_Index ci) const
{
  // ci is changed
  for (size_t d = 0; d < num_dim(); ++d)
  {
    ci[d] = coordinate_to_axis_index(p[d]);
  }
}

size_t Grid_Container::cell_index_to_linear_index( const Cell_Index ci ) const
{
  // ci is const
  int i = ci[0];
  assert( i < _axis_divisions );
  size_t offset = _axis_divisions;
  for (size_t d = 1; d < num_dim(); ++d)
  {
    assert( ci[d] < _axis_divisions );
    assert( ci[d] >= 0 );
    i += offset * ci[d];
    offset *= _axis_divisions;
  }
  assert(i < _num_cells && i >= 0);
  // consider returning a flag so the caller knows it is out of range
  return (size_t) i;
}


