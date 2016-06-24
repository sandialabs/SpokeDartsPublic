// Domain.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "Point_Tool.hpp"

class Domain
{    
public:
  // pass in points defining xmax and xmin
  Domain(size_t num_dim, bool is_periodic, const double * xmax, const double * xmin ) :
  _bounds_tol(0.), _own_xmax(0), _own_xmin(0), _xmax(xmax), _xmin(xmin), _pt(num_dim)
  { set(num_dim, is_periodic, xmax, xmin); }

  // pass in scalar defining xmax and xmin
  Domain(size_t num_dim, bool is_periodic = false, const double xmax = 1., const double xmin = 0. ) :
  _bounds_tol(0.), _own_xmax(0), _own_xmin(0), _xmax(0), _xmin(0), _pt(num_dim)
  { set(num_dim, is_periodic, xmax, xmin); }

  virtual ~Domain() { delete_ownx(); }

  // set by xmax pointer
  void set(size_t num_dim, bool is_periodic, const double * xmax, const double * xmin)
  { 
    delete_ownx();
    _num_dim = num_dim;
    _is_periodic = is_periodic;
    _pt.set_num_dim(num_dim);

    _xmax = xmax;
    _xmin = xmin;
  }

  // set by xmax scalar
  void set(size_t num_dim, bool is_periodic = false, const double xmax = 1. , const double xmin = 0.)
  {
    delete_ownx(); // be sure to delete before changing _num_dim
    _num_dim = num_dim; 
    _is_periodic = is_periodic;
    _pt.set_num_dim(num_dim);

    _own_xmax = _pt.new_point();
    _own_xmin = _pt.new_point();
    _pt.assign(_own_xmax, xmax);
    _pt.assign(_own_xmin, xmin);
    _xmax = _own_xmax;
    _xmin = _own_xmin;
  }

  void set_dimension(size_t num_dim)
  { if (num_dim != _num_dim) { set(num_dim, _is_periodic, _xmax ? _xmax[0] : 1., _xmin ? _xmin[0] : 0.); } }

  // access
  bool is_periodic() const {return _is_periodic;}
  const double * xmax() const {return _xmax;}
  const double * xmin() const {return _xmin;}
  size_t num_dim() const {return _num_dim;}

  double _bounds_tol;
  bool in_bounds_tol(const double *p) const // is p inside the non-periodic domain?
  { return _pt.in_box(p, _xmax,_xmin, _bounds_tol); }
  bool out_bounds_tol(const double *p) const // is p outside the non-periodic domain?
  { return !_pt.in_box(p, _xmax,_xmin, -_bounds_tol); }

  bool in_bounds(const double *p) const // is p inside the non-periodic domain?
  { return _pt.in_box(p, _xmax,_xmin); }
  bool in_domain(const double *p) const // is p inside the domain? always true if periodic
  {return _is_periodic || in_bounds(p); }


protected: 

  bool _is_periodic;
  double * _own_xmax, *_own_xmin;
  const double *_xmax, *_xmin;
  size_t _num_dim;
  Point_Tool _pt;

private:
  void delete_ownx()
  {
    _pt.delete_point(_own_xmax);
    _pt.delete_point(_own_xmin);
  }

};

#endif