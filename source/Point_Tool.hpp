// Point_Tool.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef POINT_TOOL_HPP
#define POINT_TOOL_HPP

#include <vector>
#include <unordered_map>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <iostream>

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif  

#ifndef E
#define E  2.7182818284590452353602874713526
#endif  


// comment the following to turn off our home-grown memory management,
// This is useful to enable debugging tools to find memory bound errors.
#define MPS_MEMORY_MANAGEMENT

class Sphere_Array
{
public:
  Sphere_Array( size_t num_dim ) : _num_dimp1(num_dim+1), _mem(0), _array_size(0) {}
  Sphere_Array( size_t num_dim, size_t array_size ) : _num_dimp1(num_dim+1), _mem(0), _array_size(array_size)
  { new_memory(_array_size); }
  ~Sphere_Array() { delete_memory(); }
  void new_memory( size_t array_size )
  {
    delete_memory();
    _array_size = array_size;
    _mem = new double[_array_size * _num_dimp1];
  }
  void delete_memory() { delete [] _mem; _mem = 0; }

  size_t size() const {return _array_size;}
  // double *sphere = sphere_array[9];
  // rvalue const access 
  const double* operator[](size_t i) const
  {
    assert(i < size());
    return & _mem[i*_num_dimp1];
  };
  // lvalue non-const access
  double * operator[](std::size_t i)
  {
    return & _mem[i*_num_dimp1];
  };

private:
  // pointer to the actual memory
  double * _mem;
  const size_t _num_dimp1;
  size_t _array_size;
};



class Point_Tool
{

public:    

  Point_Tool(size_t num_dim) : _num_dim(num_dim), _unit_disk_volume(0.) {}

  size_t num_dim() const { return _num_dim; }
  virtual void set_num_dim( size_t num_dim )
  { if (_num_dim != num_dim) { _num_dim = num_dim; _unit_disk_volume = 0.; } }
  static size_t bad_num_dim() { return size_t(-1); }
  bool valid_num_dim() const { return _num_dim != bad_num_dim(); }

  // Warning: there is nothing to prevent a "double * new_point" from being treated as a "double * new_spheres"
  // "radius" or deleting a point as a sphere will access memory allocated to some other point, or perhaps free space

  double *new_point();
  double *new_sphere(); // is a point plus radius
  Sphere_Array *new_sphere_array( size_t array_size ); // is a point plus radius
  double *new_copy_sphere( const double *s );
  double *new_copy_point( const double *p );
  void delete_point( double *&p );
  void delete_sphere( double *&s );
  void delete_sphere_array( Sphere_Array *&array ); 

  // reference to the radius of the sphere
  // uninitialized memory access if s is just a point!
  double &radius( double *s ) const { return s[_num_dim]; }
  double radius( const double *s ) const { return s[_num_dim]; }

  void multiply( double *p, const double scalar ) const;                  // *p *= scaler
  void multiply( double *p, const double *q, const double scalar ) const; // *p = scaler * q
  void add( double *p, const double *q ) const;                           // *p += *q
  void add( double *p, const double *q, const double *x ) const;          // *p = q + x
  void add( double *p, const double *q, const double a ) const;           // *p = q + scalar
  void subtract( double *p, const double *q ) const;                      // *p -= *q
  void subtract( double *p, const double *q, const double * x ) const;    // *p = q - x
  void assign_sphere( double *p, const double *q ) const;                        // *p = *q, also radius
  void assign( double *p, const double *q ) const;                        // *p = *q
  void assign( double *p, const double a ) const;                         // *p = a, each coordinate
  void axpy( double *p, const double a, const double *x,
            const double *y ) const;                                      // p = a * x + y
  void axpby( double *p, const double a, const double *x,
             const double b, const double *y ) const;                     // p = a * x + b * y
  double dot_product(const double *p, const double *q) const;             // dot product of p and q
  double norm_squared(const double *p) const;                             // p dot p
  double distance_squared(const double *p, const double *q) const;        // p-q dot p-q
  double normalize(double *p) const; // return the original norm of p
  
  // find the ghost of q that is closest to p, put it in g
  // return the distance squared to the ghost
  // is_q is true if *g == *p
  // hardcoded for a domain width of 1 in all dimensions; easy to make a version that takes box boundaries as input
  double closest_ghost( const double *p, const double *q, double *g, bool &is_q );
  double closest_ghost_distance_squared( const double *p, const double *q );
  // return the closest ghost as above, but also return the distance to the second-closest ghost
  double two_closest_ghosts( const double *p, const double *q, double *g, bool &is_q,
                            /*double &first_distance_squared,*/ double &second_distance_squared );
  

  // convert a point specified by a*u + c into coordinates p
  void ray_to_coordinates( const double *c, const double *u, double a, double *p) const
  {  axpy( p, a, u, c ); }

  // convert a point specified by a*u + b*v + c into coordinates p
  // a and b are normalized for radius 1, so neeed to specify radius as well
  void arc_to_coordinates( const double *c, const double *u, const double *v, double a, double b, double r, double *p) const
  {  axpby( p, r*a, u, r*b, v ); }

  void theta_to_xy( double theta, double &x, double &y) const
  { x = cos(theta); y = sin(theta); }

  // unit polar coordinates (angles) to x,y,z coordinates
  void theta_to_xyz( const double *theta, double *p) const;
  
  // iterate over a uniform angular grid, which is obviously not uniform on the sphere
  // start theta at all zeros
  // returns false when the last angle has been reached
  // note in d dimensions, there are only d-1 angles, so ignore the last one
  bool next_theta( double *theta, const double increment ) const;

  
  void first_quadrant_theta( double *theta );
  bool next_quadrant_theta( double *theta, const double increment ) const;



  // volume calculations
  // return the volume of a radius-r sphere, divided by the volume of a radius 1 sphere
  double relative_volume(double r) const
  {
    return pow(r, num_dim());
  }
  // for an annulus
  double relative_volume(double r_inner, double r_outer) const  
  {
    return relative_volume(r_outer) - relative_volume(r_inner);
  }
  // return the value R such that the ratio of the volume of the radius-R sphere to the radius-r sphere is t.
  double inverse_relative_volume(double t, double r) const
  {
    const double v = t * pow(r, num_dim());
    const double ir = pow(v, 1. / num_dim());
    return ir;
  }
  // return the value r such that the ratio of the volume of the annulus from r_inner to r is
  // t times the volume of the annulus form r_inner to r_outer
  double inverse_relative_volume(double t, double r_inner, double r_outer) const
  {
    const double v_inner = relative_volume(r_inner);
    const double v_outer = relative_volume(r_outer);
    const double p = t * (v_outer - v_inner);
    const double r = pow( p + v_inner, 1. / num_dim() );
    return r;
  }

  // absolute volume of a disk of dimension and radius
  double absolute_volume(const double r);
  double absolute_volume(const double *sphere)
  { return absolute_volume( radius(sphere) ); }

  double box_volume( const double *xmax, const double *xmin, double frame_size  = 0.) const;
  bool in_box(const double *p, const double *xmax, const double *xmin, double frame_size = 0.) const;

  // print in human-readable format
  // note the two forms instead of an overloaded function are for the debugger
  void print_point( const double *p ) { print_point(p, std::cout); }
  void print_point( const double *p, std::ostream &out );
  void print_sphere( const double *s ) { print_sphere(s, std::cout); }
  void print_sphere( const double *s, std::ostream &out )
  {
    print_point(s, out);
    out << " r:" << radius(s);
  }
  
  // debug
  // setting a high verification level increases the run time but does extra checks.
  //  _verification_level = 0, fastest, no checks
  //  _verification_level = 1, adds a quadratic factor to overall runtime
  //  _verification_level = 2, adds a cubic factor to overall runtime
  static size_t _verification_level;
  
protected:
  size_t _num_dim;

private:
  Point_Tool();

  double _unit_disk_volume;

//=========================================================
#ifndef MPS_MEMORY_MANAGEMENT
}; // end of class Point_Tool


// use default new and delete
inline
double *Point_Tool::new_point()
{
  assert(valid_num_dim());
  return new double[num_dim()];
}

inline
double *Point_Tool::new_sphere()
{
  assert(valid_num_dim());
  return new double[num_dim()+1]; // +1 for sphere radius
}

inline
void Point_Tool::delete_point( double *&p )
{
  delete [] p;
  p = 0;
}

inline
void Point_Tool::delete_sphere( double *&s )
{
  delete [] s;
  s = 0;
}
  
//=========================================================
// block memory and recycled memory classes
#else

  // use block memory allocate new_points and new_spheres, as we use a lot of them.
  class Mem_Block
  {
  public:
    Mem_Block( size_t mem_size )
    {
      _mem = new double[mem_size];
      _next_free = _mem;
      _last_mem = &_mem[mem_size]; // this memory location is unallocated, bad
    }
    ~Mem_Block() { delete [] _mem; }
    double *get_next_mem( size_t len )
    {
      double *current_free = _next_free;
      _next_free += len; 
      return (_next_free > _last_mem) ? 0 : current_free;
    }
    double * _mem;
    double * _next_free;
    double * _last_mem;
  };
  static size_t _block_size;
  static std::vector<Mem_Block*> _mem_blocks;
  // always succeeds
  static double *get_block_mem(size_t size);
  
  // pool of freed memory, of specific lengths
  // these are dimension specific
  typedef std::vector<double *> Point_Pool; // or a sphere pool
  typedef std::pair <size_t, Point_Pool*> Pool_Pair;
  typedef std::unordered_map<size_t, Point_Pool*> Pool_Map;
  static Pool_Map _pool_map;

  // return memory from the _point_pool, if any, or 0
  static double *get_pool_mem(size_t sz);
  // creates a pool if there isn't one
  Point_Pool *force_point_pool( size_t sz );
  
  // always succeeds; Point_Pool should call this one
  static double *get_mem(size_t sz);
  
}; // end of class Point_Tool

inline
double * Point_Tool::get_pool_mem(size_t sz)
{
  Pool_Map::iterator pmi = _pool_map.find(sz);
  if ( pmi != _pool_map.end() && pmi->second && !pmi->second->empty() )
  {
    Point_Pool *pp = pmi->second;
    double *mem = pp->back();
    pp->pop_back();
    return mem;
  }
  return 0;
}

inline
double *Point_Tool::get_mem( size_t sz )
{
  assert( sz > 0 && sz != bad_num_dim() );
  double *mem = get_pool_mem(sz);
  if (!mem)
    mem = get_block_mem(sz);
  assert(mem);
  return mem;
}

inline
Point_Tool::Point_Pool * Point_Tool::force_point_pool( size_t sz )
{
  Pool_Map::iterator pmi = _pool_map.find(sz);
  if ( pmi != _pool_map.end() )
    return pmi->second;
  Point_Pool *new_pool = new Point_Pool();
  // prefered c++11 way
  // _pool_map.emplace(sz, new_pool);
  // pre C++11 compilers don't know about emplace.
  // instead do:
  _pool_map[sz] = new_pool;
  return new_pool;
}

inline
double *Point_Tool::new_point()
{
  // return new double[num_dim()];
  assert(valid_num_dim());
  return get_mem( num_dim() );
}

inline
double *Point_Tool::new_sphere()
{
  //return new double[num_dim()+1]; // +1 for sphere radius
  assert(valid_num_dim());
  return get_mem( num_dim() + 1 );
}

inline
void Point_Tool::delete_point( double *&p )
{
  // delete [] p;
  if (p)
  {
    assert(valid_num_dim());
    Point_Pool *pp = force_point_pool( num_dim() );
    pp->push_back(p);
    p = 0;
  }
}
inline
void Point_Tool::delete_sphere( double *&s )
{
  // delete [] s;
  if (s)
  {
    assert(valid_num_dim());
    Point_Pool *pp = force_point_pool( num_dim() + 1 );
    pp->push_back(s);
    s = 0;
  }
}

#endif
// end of memory management
//=========================================================


// Sphere arrays get memory from the system
inline
Sphere_Array *Point_Tool::new_sphere_array( size_t array_size )
{
  assert(valid_num_dim());
  return new Sphere_Array( array_size, num_dim() );
}
inline
void delete_sphere_array( Sphere_Array *&array )
{
  delete array;
  array = 0;
}

inline
double *Point_Tool::new_copy_sphere(const double *s )
{
  assert(valid_num_dim());
  double *c = new_sphere();
  assign(c , s);
  radius(c) = radius(s);
  return c;
}

inline
double *Point_Tool::new_copy_point(const double *p )
{
  assert(valid_num_dim());
  double *c = new_point();
  assign(c , p);
  return c;
}

inline
void Point_Tool::assign_sphere( double *c, const double *s ) const
{
  assert(valid_num_dim());
  assign(c , s);
  radius(c) = radius(s);
}

inline
void Point_Tool::multiply( double *p, const double scalar ) const
{
  // p *= scaler
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] *= scalar;
  }
}

inline
void Point_Tool::multiply( double *p, const double *q, const double scalar ) const
{
  // *p = scaler * q
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = scalar * q[idim];
  }
}

inline
void Point_Tool::add( double *p, const double *q, const double a ) const
{
  // *p = q + a_scalar
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = q[idim] + a;
  }
}

inline
void Point_Tool::add( double *p, const double *q ) const
{
  // p += q
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] += q[idim];
  }
}

inline
void Point_Tool::add( double *p, const double *q, const double *x ) const
{
  // *p = q + x
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = q[idim] + x[idim];
  }
  
}

inline
void Point_Tool::subtract( double *p, const double *q ) const
{
  // *p -= *q
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] -= q[idim];
  }
}

inline
void Point_Tool::subtract( double *p, const double *q, const double * x ) const
{
  // *p = q - x
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = q[idim] - x[idim];
  }
}

inline
void Point_Tool::assign( double *p, const double *q ) const
{
  // *p = *q
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = q[idim];
  }
}

inline
void Point_Tool::assign( double *p, const double a ) const
{
  // *p = a, each coordinate
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = a;
  }
}

inline
void Point_Tool::axpy( double *p, const double a, const double *x, const double *y ) const
{
  // p = a * x + y
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = a * x[idim] + y[idim];
  }
}

inline
void Point_Tool::axpby( double *p, const double a, const double *x,
                               const double b, const double *y ) const
{
  // p = a * x + b * y
  assert(valid_num_dim());
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    p[idim] = a * x[idim] + b * y[idim];
  }
}

inline
double Point_Tool::dot_product(const double *p, const double *q) const
{
  // dot product of p and q
  assert(valid_num_dim());
  double dot(0.);
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    dot += p[idim] * q[idim];
  }
  return dot;
}

inline
double Point_Tool::norm_squared(const double *p) const
{
  // p dot p
  assert(valid_num_dim());
  double n(0.);
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    n += p[idim] * p[idim];
  }
  return n;
}

inline
double Point_Tool::distance_squared(const double *p, const double *q) const
{
  // p dot p
  assert(valid_num_dim());
  double d(0.);
  for ( size_t idim = 0; idim < _num_dim; ++idim)
  {
    double di(p[idim] - q[idim]);
    d += di*di;
  }
  return d;
}

inline
double Point_Tool::normalize(double *p) const
{
  const double norm = sqrt(norm_squared(p));
  multiply(p, 1. / norm );
  return norm;
}


inline
bool Point_Tool::in_box(const double *p, const double *xmax, const double *xmin, double frame_size) const
{
  assert(valid_num_dim());
  for (size_t d = 0; d < num_dim(); ++d)
  {
    if (p[d] < xmin[d] - frame_size || p[d] >= xmax[d] + frame_size)
      return false;
  }
  return true;
}

inline
double Point_Tool::closest_ghost( const double *p, const double *q, double *g, bool & is_q)
{
  is_q = true;
  double dx2 = 0.; // sum of squared coordinate differences
  radius(g) = radius(q);
  for (size_t d = 0; d < num_dim(); ++d)
  {
    double dx = p[d]-q[d];
    g[d] = q[d];
    if ( dx > 0.5 ) // 1/2 domain width, hardcoded
    {
      do
      {
        g[d]++;
        dx--;
      } while ( dx > 0.5 );
      is_q = false;
    }
    else if ( dx < -0.5 )
    {
      do
      {
        g[d]--;
        dx++;
      } while ( dx < -0.5 );
      is_q = false;
    }
    
    dx2 += dx * dx;
  }
  return dx2;
}

inline
double Point_Tool::closest_ghost_distance_squared( const double *p, const double *q )
{
  double dx2 = 0.; // sum of squared coordinate differences
  for (size_t d = 0; d < num_dim(); ++d)
  {
    double dx = p[d]-q[d];
    dx -= round(dx); // assumes domain width = 1
    assert( dx < 0.501 );
    assert( dx > -0.501 );
    dx2 += dx * dx;
  }
  return dx2;
}

inline
double Point_Tool::two_closest_ghosts( const double *p, const double *q, double *g, bool &is_q,
                                      double &second_distance_squared )
{
  is_q = true;
  double dx2 = 0.; // sum of squared coordinate differences
  radius(g) = radius(q);
  
  // for second_distance
  double best_dx( 10. ), replace_dx(0.); // bigger than period
  
  for (size_t d = 0; d < num_dim(); ++d)
  {
    double dx = p[d]-q[d];
    g[d] = q[d];
    if ( dx > 0.5 ) // 1/2 domain width, hardcoded
    {
      do
      {
        g[d]++;
        dx--;
      } while ( dx > 0.5 );
      is_q = false;
    }
    else if ( dx < -0.5 )
    {
      do
      {
        g[d]--;
        dx++;
      } while ( dx < -0.5 );
      is_q = false;
    }

    dx2 += dx * dx;
    
    // how does the distance change if we used the second-closest periodic copy of this coordinate?
    const double second_dx = 1. - fabs(dx);
    // if this is the smallest change so far, save it
    if (second_dx < best_dx )
    {
      best_dx = second_dx;
      replace_dx = dx;
    }
  }
  second_distance_squared = dx2 + best_dx * best_dx - replace_dx * replace_dx;
  
  return dx2;
}

inline
void Point_Tool::theta_to_xyz( const double *theta, double *p) const
{
  p[0] = cos(theta[0]);
  for (size_t d = 1; d < num_dim(); ++d)
  {
    double v = 1.;
    for (size_t dd = 0; dd < d; ++dd )
    {
      v *= sin(theta[dd]);
    }
    if ( d < num_dim()-1 )
      p[d] = v * cos(theta[d]);
    else
      p[d] = v * sin(theta[d-1]);
  }
}


inline
bool Point_Tool::next_theta( double *theta, const double increment ) const
{
  //theta[0] in 0, 2pi
  //theta[...] in 0, pi
  theta[0] += increment;
  if (theta[0] > 2*PI)
  {
    theta[0] = 0.;
    theta[1] += increment;
    for (size_t d = 1; d < num_dim()-1; ++d)
    {
      if (theta[d] > PI)
      {
        theta[d] = 0.;
        theta[d+1] += increment;
      }
      else
        break;
    }
  }
  
//  // debug
//  std::cout << " theta = ";
//  for (size_t dd = 0; dd < num_dim(); ++dd )
//    std::cout << theta[dd] << " ";
//  std::cout << std::endl;
  
  return (theta[num_dim()-1] == 0.);
}


inline
void Point_Tool::first_quadrant_theta( double *theta )
{
  assign(theta, PI/4.);
  theta[num_dim()-1] = 0.;
}

inline
bool Point_Tool::next_quadrant_theta( double *theta, const double increment ) const
{
  //theta[0] in 0, 2pi:  go over pi/4 to 3pi/4
  //theta[...] in 0, pi:  go over pi/4 to 3pi/4
  theta[0] += increment;
  if (theta[0] > 3.*PI/4.)
  {
    theta[0] = PI/4.;
    theta[1] += increment;
    for (size_t d = 1; d < num_dim()-1; ++d)
    {
      if (theta[d] > 3.*PI/4.)
      {
        theta[d] = PI/4.;
        theta[d+1] += increment;
      }
      else
        break;
    }
  }
  
  //  // debug
  //  std::cout << " theta = ";
  //  for (size_t dd = 0; dd < num_dim(); ++dd )
  //    std::cout << theta[dd] << " ";
  //  std::cout << std::endl;
  
  return (theta[num_dim()-1] == 0.);
}



#endif