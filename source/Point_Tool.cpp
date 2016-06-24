// Point_Tool.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Point_Tool.hpp"

size_t Point_Tool::_verification_level(0);

//=============================================
#ifdef MPS_MEMORY_MANAGEMENT

size_t Point_Tool::_block_size( 65536 );
std::vector<Point_Tool::Mem_Block*> Point_Tool::_mem_blocks;

Point_Tool::Pool_Map Point_Tool::_pool_map;

double *Point_Tool::get_block_mem(size_t size)
{
  double *mem(0);
  if (!_mem_blocks.empty())
    mem = _mem_blocks.back()->get_next_mem(size);
  if (!mem)
  {
    _block_size *= 2; // or 1.5?
    if (_block_size < size)
      _block_size = size;
    _mem_blocks.push_back( new Mem_Block( _block_size ) );
    mem = _mem_blocks.back()->get_next_mem(size);
  }
  assert(mem);
  return mem;
}
#endif
//=============================================

double Point_Tool::absolute_volume(const double r)
{
  // use saved unit volume from last time
  if ( _unit_disk_volume == 0. )
  {
    // "0-dim"
    double v = 1.; // volume of unit ball
    double s = 2.; // surface area of one dim higher ball
    // to reduce numerical issues, huge number divided by even huger number,
    // we use iteration
    for (unsigned int i = 1; i <= num_dim(); ++i)
    {
      double v_next = s / i;
      s = 2 * PI * v;
      v = v_next;
    }
    _unit_disk_volume = v;
  }
  
  const double volume = _unit_disk_volume * relative_volume(r);
  
  return volume;
}


double Point_Tool::box_volume( const double *xmax, const double *xmin, double frame_size ) const
{
  double vol(1.);
  for (size_t d = 0; d < num_dim(); ++d)
    vol *= xmax[d] - xmin[d] + 2 * frame_size;
  return vol;
}

void Point_Tool::print_point( const double *p, std::ostream &out )
{
  if (!p)
  {
    out << "NULL";
    return;
  }
  out << "[ " << p[0];
  for (size_t d = 1; d < num_dim(); ++d)
  {
    out << ", " << p[d];
  }
  out << " ]";
}
