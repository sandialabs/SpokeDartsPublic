// Search_Range.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// range tree
#include "Search_Range.hpp"

inline
void Search_Range::old_lo_hi( size_t d, double& x_lo, double &x_hi )
{
  Near_Workspace &ws = near_workspace;
  x_lo = ws._p[d] - ws._thresh;
  x_hi = ws._p[d] + ws._thresh;
}

bool Search_Range::new_lo_hi( size_t &node, size_t d, double &x_lo, double &x_hi)
{
  Near_Workspace &ws = near_workspace;
  
  if (ws._p == _spheres[node])
    return false;
  
  const double pair_distance = _tt.distance_squared( ws._p, _spheres[node]);

  // new nearest?
  if (pair_distance <= ws._thresh_squared)
  {
    ws._thresh_squared = pair_distance;
    ws._thresh = sqrt( pair_distance );
    ws._conflict_sphere = node;
    old_lo_hi(d, x_lo, x_hi);
    return true;
  }  
  return false;
}

inline
bool Search_Range::new_old_lo_hi( size_t &node, size_t d, double &x_lo, double &x_hi)
{
  if (new_lo_hi(node, d, x_lo, x_hi))
    return true;
  old_lo_hi(d, x_lo, x_hi);
  return false;
}

double Search_Range::nearest_sphere( const double *p, size_t &near_sphere, const double distance_threshold )
{
  _stats.start_clock();

  Near_Workspace &ws = near_workspace;
  ws._p = p;
  if ( distance_threshold < sqrt( std::numeric_limits<double>::max() - 10. ) )
  {
    ws._thresh = distance_threshold;
    ws._thresh_squared = ws._thresh * ws._thresh;
  }
  else
  {
    ws._thresh =   sqrt( std::numeric_limits<double>::max() - 10. );
    ws._thresh_squared = std::numeric_limits<double>::max() - 10.;
  }
  ws._sub_container = 0;

  ws._sub_container = 0;
  ws._num_calls = 0;

  ws._conflict_sphere = Spheres::bad_sphere_index(); // stand in for searching from p

  nearest_sphere_recurse( &_tree );

  near_sphere = ws._conflict_sphere;

  _stats.collect_stats(size(), ws._num_calls, 1);

  return ws._thresh;
}
void Search_Range::nearest_sphere_recurse( One_Tree *tree )
{
  Near_Workspace &ws = near_workspace;

  // debug
  const bool debug = false;
  assert(ws._num_calls <  num_nodes() );
  ++ws._num_calls;
  
  const size_t d = tree->_d;
  double x_lo, x_hi;
  old_lo_hi(d, x_lo, x_hi);

  // walk until the path to x_lo and x_hi diverge
  // this might update x_lo and x_hi with a new closest
  Tree_Node *n = find_split( tree, x_lo, x_hi, true );

  if (debug)
  {
    std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
    std::cout << "Diverges at node ";
    Range_Tree::print_single_node( n, tree );
    // std::cout << std::endl;
  }

  // never split, nothing is in range
  if (!n)
  {
    if (debug)
      std::cout << " no nearest sphere was in range." << std::endl;
    return;
  }

  // If n is a leaf, then immediately recurse by dimension on its subtree.
  // For efficiency we just check the distance immediately, since the result is the same.
  if (n->is_leaf())
  {
    if (debug)
      std::cout << " fathomed at leaf node." << std::endl;
    new_lo_hi( n->_sphere, d, x_lo, x_hi);
    return;
  }

  // go down right and left sides, checking subtrees
  // right
  {
    if (debug)
    {
      std::cout << "  Searching Right subtree " << std::endl;
    }
    Tree_Node * r = n->_right;
    while( r )
    {
      new_old_lo_hi( r->_sphere, d, x_lo, x_hi);

      const double &x = _spheres[ r->_sphere ][ d ];
      
      if ( r->is_leaf() )
        break;

      // discard right->right?
      if ( x > x_hi )
        r = r->_left;
      else
      {
        // recurse on left
        if (r->_left->_subtree)
        {
          if (debug) std::cout << "recursively searching left subtree. " << std::endl;
          nearest_sphere_recurse( r->_left->_subtree );
        }
        // fathomed at last dimension
        else
        {
          if (debug)
          {
            std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
            std::cout << " Nearest all leaves, left subtree, rooted at node ";
            print_single_node(r->_left, tree);
            std::cout << std::endl;
          }
          assert( last_dimension(tree) || r->_left->is_leaf() );
          // go down left (recurse) but with current dimension
          nearest_sphere_leaves( r->_left, d );
        }
        
        r = r->_right;
      }
    }
  }
  // left
  {
    if (debug)
    {
      std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
      std::cout << "Searching Left subtree " << std::endl;
    }
    Tree_Node * lft = n->_left;
    while( lft )
    {
      new_old_lo_hi( lft->_sphere, d, x_lo, x_hi);
      
      if (lft->is_leaf())
        break;

      const double &x = _spheres[ lft->_sphere ][ tree->_d ];
      if ( x_lo > x )
        lft = lft->_right;
      else
      {
        // right tree is (was) in the range
        if ( lft->_right )
        {
          if ( lft->_right->_subtree )
            nearest_sphere_recurse( lft->_right->_subtree );
          // fathomed at last dimension
          else
          {
            if (debug)
            {
              std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
              std::cout << " Nearest all leaves, right subtree, rooted at node ";
              print_single_node(lft->_right, tree);
              std::cout << std::endl;
            }
            assert( last_dimension(tree) || lft->_right->is_leaf() );
            // go down right (recurse) but with current dimension
            nearest_sphere_leaves( lft->_right, d );
          }
        }
        lft = lft->_left;
        // null if it was a leaf
      }
    }
  }
}


void Search_Range::nearest_sphere_leaves( Tree_Node *node, size_t d )
{
  double x_lo, x_hi;
  new_old_lo_hi(node->_sphere, d, x_lo, x_hi);

  if (node->is_leaf())
    return;
  
  Tree_Node *lft = node->_left;
  Tree_Node *rgh = node->_right;
  new_lo_hi(lft->_sphere, d, x_lo, x_hi);
  new_lo_hi(rgh->_sphere, d, x_lo, x_hi);
  
  // test, exhaustive recursion
//  nearest_sphere_leaves(lft, d );
//  nearest_sphere_leaves(rgh, d);
//  return;

  const double x = _spheres[node->_sphere][d];
  const double xl = _spheres[node->_left->_sphere][d];
  const double xr = _spheres[node->_right->_sphere][d];

  if (!lft->is_leaf())
  {
    if (x >= x_lo && xl <= x_hi)
    {
      nearest_sphere_leaves(lft->_right, d);
      old_lo_hi(d, x_lo, x_hi);
    }
    if (xl >= x_lo)
    {
      nearest_sphere_leaves(lft->_left, d);
      old_lo_hi(d, x_lo, x_hi);
    }
  }
  if (!rgh->is_leaf())
  {
    if (xr >= x_lo && x <= x_hi)
    {
      nearest_sphere_leaves(rgh->_left, d);
      old_lo_hi(d, x_lo, x_hi);
    }
    if (xr <= x_hi)
    {
      nearest_sphere_leaves(rgh->_right, d);
      old_lo_hi(d, x_lo, x_hi);
    }
  }
}


void Search_Range::all_near_spheres(Sphere_Container *sub_container, const double *p, const double dist, const bool subtract_radius )
{
  _stats.start_clock();
  sub_container->clear();
  if (!size())
    return;

  // stuff that would normally be passed to the recursive function
  Near_Workspace &ws = near_workspace;
  ws._p = p;
  ws._thresh = dist;
  if (subtract_radius)
    ws._thresh += _spheres.radius(0);
  ws._thresh_squared = ws._thresh * ws._thresh;

  ws._sub_container = sub_container;
  ws._num_calls = 0;

  
  // this function strongly affects the overall run-time
  all_near_spheres_recurse( &_tree );

  _stats.collect_stats(size(), ws._num_calls, sub_container->size());
}
void Search_Range::all_near_spheres_recurse( One_Tree *tree )
{
  const bool debug = false;
  
  Near_Workspace &ws = near_workspace;

  // debug
  assert(ws._num_calls <  num_nodes() );
  ++ws._num_calls;
  
  const size_t d = tree->_d;
  double x_lo = ws._p[d] - ws._thresh;
  double x_hi = ws._p[d] + ws._thresh;

  // walk until the path to x_lo and x_hi diverge
  Tree_Node *n = find_split( tree, x_lo, x_hi );
  if (debug)
  {
    std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
    std::cout << "Diverges at node ";
    Range_Tree::print_single_node( n, tree );
    // std::cout << std::endl;
  }

  // never split, nothing is in range
  if (!n)
  {
    if (debug)
      std::cout << " no spheres were in range." << std::endl;
    return;
  }

  // If n is a leaf, then immediately recurse by dimension on its subtree.
  // For efficiency we just check the distance immediately, since the result is the same.
  if (n->is_leaf())
  {
    if (debug)
      std::cout << " fathomed at leaf node." << std::endl;
    add_leaf(n->_sphere);
    return;
  }

    // go down right and left sides, reporting subtrees
  // right
  {
    if (debug)
    {
      std::cout << "  Searching Right subtree " << std::endl;
    }
    Tree_Node * r = n->_right;
    while( r )
    {
      const double &x = _spheres[ r->_sphere ][ d ];
      if ( x_hi > x )
      {
        // left tree is in the range
        if ( r->_left )
        {
          if (r->_left->_subtree)
            all_near_spheres_recurse( r->_left->_subtree );
          // fathomed at last dimension
          else
          {
            if (debug)
            {
              std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
              std::cout << " Adding left subtree rooted at node ";
              print_single_node(r->_left, tree);
              std::cout << std::endl;
            }
            assert( last_dimension(tree) || r->_left->is_leaf() );
            add_all_leaves( r->_left );
          }
        }
        // leaf
        else
        {
          assert( r->is_leaf() );
          if (debug)
          {
            std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
            std::cout << " Adding leaf node ";
            print_single_node(r, tree);
            std::cout << std::endl;
          }
          add_leaf( r->_sphere );
        }
        r = r->_right;
        // null if it was a leaf
      }
      else
        r = r->_left;
        // if r == 0, then its parent was a leaf and out of range
    }
  }
  // left
  {
    if (debug)
    {
      std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
      std::cout << "Searching Left subtree " << std::endl;
    }
    Tree_Node * lft = n->_left;
    while( lft )
    {
      const double &x = _spheres[ lft->_sphere ][ tree->_d ];
      if ( x_lo <= x )
      {
        // right tree is in the range
        if ( lft->_right )
        {
          if ( lft->_right->_subtree )
            all_near_spheres_recurse( lft->_right->_subtree );
          // fathomed at last dimension
          else
          {
            if (debug)
            {
              std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
              std::cout << " Adding right subtree rooted at node ";
              print_single_node(lft->_right, tree);
              std::cout << std::endl;
            }
            assert( last_dimension(tree) || lft->_right->is_leaf() );
            add_all_leaves( lft->_right );
          }
        }
        // leaf
        else
        {
          assert( lft->is_leaf() );
          if (debug)
          {
            std::cout << "Search_Range [" << x_lo << ", " << x_hi << "] coordinate " << tree->_d << ". ";
            std::cout << " Adding leaf node ";
            print_single_node(lft, tree);
            std::cout << std::endl;
          }
          add_leaf( lft->_sphere );
        }
        lft = lft->_left; 
        // null if it was a leaf
      }
      else
      {
        lft = lft->_right;
        // if lft == 0, then its old value (parent) was a leaf and out of range
      }
    }
  }
}

void Search_Range::add_all_leaves( Tree_Node *node )
{
  // for efficiency, this should only be called after fathoming by dimension
  if ( node->is_leaf() )
  {
    add_leaf( node->_sphere );
  }
  else
  {
    add_all_leaves( node->_left );
    add_all_leaves( node->_right );
  }
}

void Search_Range::add_leaf( size_t isphere )
{
  assert( _spheres.is_valid_sphere(isphere));
  const double *s = _spheres[isphere];
  Near_Workspace &ws = near_workspace;
  if (s != ws._p)
  {
    const double pair_distance_squared = _tt.distance_squared( ws._p, s );
    // if close add to list of found vertices
    if (pair_distance_squared < ws._thresh_squared)
    {
      ws._sub_container->add_sphere(isphere);
    }
  }
}




bool Search_Range::no_near_spheres( const double *p, size_t &near_sphere, const double radius_factor, const double dist)
{
  _stats.start_clock();
  if (!size())
  {
    near_sphere = _spheres.bad_sphere_index();
    return true;
  }

  // stuff that would normally be passed to the recursive function
  Near_Workspace &ws = near_workspace;
  ws._p = p;
  ws._thresh = radius_factor * _spheres.radius(0)  + dist - 1e-5;
  ws._thresh_squared = ws._thresh * ws._thresh;

  ws._sub_container = 0;
  ws._num_calls = 0;

  ws._conflict_sphere = _spheres.bad_sphere_index();

  bool return_value = no_near_spheres_recurse( &_tree );
  near_sphere = ws._conflict_sphere;
  
  _stats.collect_stats(size(), ws._num_calls, !return_value);

  return return_value;
}
bool Search_Range::no_near_spheres_recurse( One_Tree *tree ) 
{

  Near_Workspace &ws = near_workspace;

  // debug
  assert(ws._num_calls <= num_nodes() );
  ++ws._num_calls;
  
  double x_lo = ws._p[tree->_d] - ws._thresh;
  double x_hi = ws._p[tree->_d] + ws._thresh;

  // walk until the path to x_lo and x_hi diverge
  Tree_Node *n = find_split( tree, x_lo, x_hi);

  // nothing in range
  if (!n)
    return true;

  // only one sphere in range, check it now
  if (n->is_leaf())
    return no_near_leaf(n->_sphere);

  // go down right and left sides
  // right
  {
    Tree_Node * r = n->_right;
    while( r )
    {
      assert( r->is_wellformed() );
      const double &x = _spheres[ r->_sphere ][ tree->_d ];
      if ( x_hi > x )
      {
        // left tree is in the range
        if ( r->_left )
        {
          assert( r->_left->is_wellformed() );
          if (r->_left->_subtree)
          {
            if (!no_near_spheres_recurse( r->_left->_subtree ))
              return false;
          }
          // fathomed at last dimension
          else
          {
            assert( last_dimension(tree)  );
            // check entire subtree explicitly
            if (!no_near_leaves( r->_left ))
              return false;
          }
        }
        // leaf
        else
        {
          // check leaf sphere
          assert(r->is_leaf());
          if (!no_near_leaf( r->_sphere ))
            return false;
        }
        r = r->_right;
        // null if it was a leaf
      }
      else
        r = r->_left;
        // if r == 0, then its parent was a leaf and out of range
    }
  }
  // left
  {
    Tree_Node * lft = n->_left;
    while( lft )
    {
      assert( lft->is_wellformed() );
      const double &x = _spheres[ lft->_sphere ][ tree->_d ];
      if ( x_lo <= x )
      {
        // right tree is in the range
        if ( lft->_right )
        {
          assert( lft->is_wellformed() );
          if ( lft->_right->_subtree )
          {
            if (!no_near_spheres_recurse( lft->_right->_subtree ))
              return false;
          }
          // fathomed at last dimension
          else
          {
            assert( last_dimension(tree)  );
            if (!no_near_leaves( lft->_right ))
              return false;
          }
        }
        // leaf
        else
        {
          assert( lft->is_leaf() );
          if (!no_near_leaf( lft->_sphere ))
            return false;
        }
        lft = lft->_left; 
        // null if it was a leaf
      }
      else
      {
        lft = lft->_right;
        // if lft == 0, then its old value (parent) was a leaf and out of range
      }
    }
  }
  return true;
}

bool Search_Range::no_near_leaves( Tree_Node *node )
{
  if ( node->is_leaf() )
    return no_near_leaf( node->_sphere );

  // non-leaf
  if (!no_near_leaves( node->_left ) || !no_near_leaves( node->_right ))
    return false;
  return true;
}

bool Search_Range::no_near_leaf( size_t isphere )
{
  assert( _spheres.is_valid_sphere(isphere) );
  const double *s = _spheres[isphere];
  Near_Workspace &ws = near_workspace;
  if (ws._p != s)
  {
    const double pair_distance = _tt.distance_squared( ws._p, s);

    if (pair_distance < ws._thresh_squared) // on surface should be OK
    {
      ws._conflict_sphere = isphere;
      return false; 
    }
  }
  return true;
}



void Search_Range::trim_line_anchored(const double* c, double* u, double r, double A, double &A_1, double &A_2, const double radius_factor )
{
  _stats.start_clock();
  if (!size())
    return;
  Trim_Statespace * z = &trim_statespace;
  z->set(c, u, r, A, A_1, A_2, radius_factor);
  z->_visited = 0; z->_found = 0;
  trim_line_anchored_recurse( &_tree );
  _stats.collect_stats(size(), z->_visited, z->_found);
}

void Search_Range::get_line_neighborhood( size_t d, double &x_lo, double &x_hi )
{
  Trim_Statespace * z = &trim_statespace;
  {
    const double x_1 = z->_c[d] + *(z->_A_1) * z->_u[d];   // "x" coord of point at A_1
    const double x_2 = z->_c[d] + *(z->_A_2) * z->_u[d];   // "x" coord of point at A_2
    if ( x_1 < x_2 )
    {
      x_lo = x_1; 
      x_hi = x_2;
    }
    else
    {
      x_lo = x_2;
      x_hi = x_1;
    }
  }
  x_lo -= z->_r;
  x_hi += z->_r;
}

void Search_Range::trim_line_anchored_recurse( One_Tree *tree )
{

  // walk until the path to x_lo and x_hi diverge
  double x_lo, x_hi;
  Tree_Node *n = find_split( tree, x_lo, x_hi );

  // nothing was in range, no spheres near enough to possibly trim the line segment
  if (!n)
    return;

  // leaf, only one sphere in range, trim with it now
  if (n->is_leaf())
  {
    trim_leaf( n->_sphere );
    return;
  }

   // go down right and left sides
  // right
  {
    Tree_Node * r = n->_right;
    while( r )
    {
      assert( r->is_wellformed() );
      const double &x = _spheres[ r->_sphere ][ tree->_d ];
      if ( x_hi > x )
      {
        // left tree is in the range
        if ( r->_left )
        {
          assert( r->_left->is_wellformed() );
          if (r->_left->_subtree)
            trim_line_anchored_recurse( r->_left->_subtree );
          // fathomed at last dimension
          else
          {
            assert( last_dimension(tree)  );
            trim_all_leaves( r->_left );
            get_line_neighborhood( tree->_d, x_lo, x_hi );
          }
        }
        // leaf
        else
        {
          assert( r->is_leaf() );
          trim_leaf( r->_sphere );
          get_line_neighborhood( tree->_d, x_lo, x_hi );
        }
        r = r->_right;
        // null if it was a leaf
      }
      else
        r = r->_left;
        // if r == 0, then its parent was a leaf and out of range
    }
  }
  // left
  {
    Tree_Node * lft = n->_left;
    while( lft )
    {
      assert( lft->is_wellformed() );
      const double &x = _spheres[ lft->_sphere ][ tree->_d ];
      if ( x_lo <= x )
      {
        // right tree is in the range
        if ( lft->_right )
        {
          assert( lft->_right->is_wellformed() );
          if ( lft->_right->_subtree )
            trim_line_anchored_recurse( lft->_right->_subtree );
          // fathomed at last dimension
          else
          {
            assert( last_dimension(tree)  );
            trim_all_leaves( lft->_right );
            get_line_neighborhood( tree->_d, x_lo, x_hi );
          }
        }
        // leaf
        else
        {
          assert( lft->is_leaf() );
          trim_leaf( lft->_sphere );
          get_line_neighborhood( tree->_d, x_lo, x_hi );
        }
        lft = lft->_left; 
        // null if it was a leaf
      }
      else
      {
        lft = lft->_right;
        // if lft == 0, then its old value (parent) was a leaf and out of range
      }
    }
  }
}

void Search_Range::trim_all_leaves( Tree_Node *node )
{
  assert( node->is_wellformed());
  if ( node->is_leaf() )
    trim_leaf( node->_sphere );
  else
  {
    trim_all_leaves( node->_left );
    trim_all_leaves( node->_right );
  }
}

void Search_Range::trim_leaf( size_t isphere )
{
  assert( _spheres.is_valid_sphere(isphere) );
  // the geometry work
  trim_statespace.anchored_trim(_tt, _spheres[isphere]);
}            


Search_Range::Tree_Node * Search_Range::find_split( One_Tree *tree, double &x_lo, double &x_hi, bool update_closest )
{
  Tree_Node *n = tree->_root;
  while (1)
  { 
    assert( n->is_wellformed() );
 
    if (update_closest)
      new_lo_hi( n->_sphere, tree->_d, x_lo, x_hi );

    const double &x = _spheres[ n->_sphere ][ tree->_d ];
    if ( (x_lo <= x && x_hi > x ) ) 
    {
      // split starts here
      return n;
    }
    // else both on same side

    // nothing is in range, node was never split
    if ( n->is_leaf() ) 
      return 0;

    // walk left or right
    if ( x_lo <= x )
      n = n->_left;
    else
      n = n->_right;
  }
}
 
