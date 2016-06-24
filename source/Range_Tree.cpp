// Range_Tree.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

#include "Range_Tree.hpp"
#include <assert.h>

#ifdef MPS_MEMORY_MANAGEMENT
Mem_Pool< Range_Tree::Tree_Node >   Range_Tree::Tree_Node::_mem_pool( 10000 ); // block size
Mem_Pool< Range_Tree::One_Tree >    Range_Tree::One_Tree::_mem_pool ( 1000 ); // smaller block size
#endif

#ifdef DEBUG_RANGE_TREE
size_t Range_Tree::One_Tree::_max_id(0);
#endif

void Range_Tree::new_memory()
{
  // delete old memory, if any
  delete_memory();
  _tree_ispheres_size = _array.max_size();
  _tree_ispheres = new size_t* [num_dim()];
  for (size_t d = 0; d < num_dim(); ++d)
    _tree_ispheres[d] = new size_t[_tree_ispheres_size];
}

void Range_Tree::delete_memory()
{
  clear_tree();
  if (_tree_ispheres)
  {
    for (size_t d = 0; d < num_dim(); ++d)
      delete [] _tree_ispheres[d];
    delete [] _tree_ispheres;
    _tree_ispheres = 0;
  }
  _tree_ispheres_size = 0;
}

void Range_Tree::clear()

{
  clear_tree();
  _array.clear();
}
void Range_Tree::clear_tree()
{
  if (_tree._root)
  {
    delete_tree_node( _tree._root );
    _tree._root = 0;
#ifdef DEBUG_RANGE_TREE
    _tree._max_node_id = 0;
#endif
  }
}

Range_Tree::One_Tree * Range_Tree::new_one_tree( size_t d )
{
  assert( d < num_dim() );
  return new One_Tree(d);
}
void Range_Tree::delete_one_tree( One_Tree *tree )
{
  if (tree->_root)
    delete_tree_node( tree->_root );
  delete tree;
  // caller should set tree = 0, if desired
}

Range_Tree::Tree_Node * Range_Tree::new_tree_node(One_Tree * tree)
{
  Tree_Node * node = new Tree_Node;
#ifdef DEBUG_RANGE_TREE
  node->_id = ++tree->_max_node_id;
  node->_level = tree->_current_level;
#endif
  return node;
};
void Range_Tree::delete_tree_node( Tree_Node *node )
{
  if (node->_left)
    delete_tree_node(node->_left);
  if (node->_right)
    delete_tree_node(node->_right);
  if (node->_subtree)
    delete_one_tree( node->_subtree );
  delete node;
  // caller should set node = 0, if desired
};

void Range_Tree::rebalance()
{
  if ( !size() )
    return;
  
  // debug
  const double past_time = _rebalance_stats.cpu_time;
  
  _rebalance_stats.start_clock();
  _rebalance_stats.num_rebalances++;
  
  clear_tree(); // leave array intact, it is used below
  if (_array.size())
    rebalance_tree( &_tree, &_array[0], _array.size());

  _rebalance_stats.collect_stats();
  const double elapsed_time = _rebalance_stats.cpu_time - past_time;
  
  // std::cout << "Rebalance Range_Tree " << _rebalance_stats.num_rebalances << " took " << elapsed_time << " seconds." << std::endl;
}

struct xd_less
{
  xd_less(size_t d, Spheres &spheres) : _d(d), _spheres(spheres) {}
  size_t _d;
  Spheres &_spheres;
  inline bool operator() (const size_t &pi, const size_t &qj)
  {
    return ( _spheres[pi][_d] <  _spheres[qj][_d] );
  }
};

size_t xd_less_d(0);
Spheres *xd_less_spheres(0);
int xd_less_fn(const void *pi, const void *qj)
{
//  const size_t i = * reinterpret_cast<const size_t *>(pi);
//  const size_t j = * reinterpret_cast<const size_t *>(qj);
//  spheres[i][d] - spheres[j][d]
  const double coord_diff = (*xd_less_spheres)[ * reinterpret_cast<const size_t *>(pi) ][xd_less_d] - (*xd_less_spheres)[ * reinterpret_cast<const size_t *>(qj) ][xd_less_d];
  
//  if ( coord_diff < 0 )
//    return -1;
//  if ( coord_diff > 0)
//    return 1;
//  return 0;
  return (coord_diff < 0.) ? -1 : ( coord_diff > 0. ? 1 : 0);
}

void Range_Tree::rebalance_tree(One_Tree *tree, const size_t *ispheres, size_t isize)
{
  // the data are sorted by some other coordinate
  // copy and sort by this tree's coordinate

  // copy
  //size_t *tree_ispheres = new size_t[isize];
  size_t *tree_ispheres = _tree_ispheres[ tree->_d ]; // reused
  assert( isize < _tree_ispheres_size );
  // it is OK to overright the dth array, since we recurse by dimension and the calling trees were for a strictly smaller dimension
  const size_t *hip = ispheres + isize;
  std::copy( ispheres, hip, tree_ispheres );

  // sort on the d_th coordinate
  // both of these implementations work, but the qsort one is faster
  
//  xd_less tree_less(tree->_d, _spheres);
//  size_t *tree_top = tree_ispheres + isize;
//  std::sort( tree_ispheres, tree_top, tree_less );

  xd_less_d = tree->_d;
  xd_less_spheres = &_spheres;
  std::qsort( tree_ispheres, isize, sizeof(size_t), xd_less_fn );

#ifdef DEBUG_RANGE_TREE
  tree->_current_level = 0;
#endif
  
  rebalance_node( tree, tree->_root, tree_ispheres, isize);
  // delete [] tree_ispheres;
}


// add ispheres[0] to [hi] to the node
// these are already sorted by the _d coordinate
void Range_Tree::rebalance_node(One_Tree *tree, Tree_Node *&node, size_t *ispheres, size_t isize)
{

  // make a new node
  assert( !node );
  node = new_tree_node(tree);

#ifdef DEBUG_RANGE_TREE
  node->_level = tree->_current_level;
#endif

  // find median isphere, assign it to node
  assert( size() );
  const size_t mid = (isize-1) / 2;
  node->_sphere = ispheres[mid];

  // recurse right and left
  // if size==1, then the node is a leaf
  if (isize > 1)
  {
#ifdef DEBUG_RANGE_TREE
    tree->_current_level++;
#endif
    rebalance_node( tree, node->_left,   ispheres,        mid+1 );
    rebalance_node( tree, node->_right, &ispheres[mid+1], isize - (mid+1) );
#ifdef DEBUG_RANGE_TREE
    tree->_current_level--;
#endif
  }

  // recurse by dimension
  assert( !node->_subtree );
  if ( tree->_d + 1 < num_dim() && isize > 1)
  {
    node->_subtree = new_one_tree( tree->_d + 1 );
    rebalance_tree( node->_subtree, ispheres, isize );
  }

}

void Range_Tree::add_spheres(Sphere_Container *spheres)
{
  // add to array
  for (size_t i = spheres->first(); i != spheres->bad_sphere_index(); i = spheres->next())
    _array.add_sphere(i);
  
  rebalance();
}

void Range_Tree::add_sphere(const size_t isphere)
{  
  // add to array
  _array.add_sphere(isphere);

  add_sphere_to_tree( isphere );
}

void Range_Tree::add_sphere_to_tree(const size_t isphere)
{
  // add to tree
  assert( _spheres.is_valid_sphere( isphere ) );
  add_sphere_to_one_tree( &_tree, isphere );
}


void Range_Tree::add_sphere_to_one_tree(One_Tree *tree, const size_t isphere)
{
  tree->_traversed_depth = 0;
  tree->_current_level = 0;
  add_sphere_to_node( tree, tree->_root, isphere );
  if (tree->_traversed_depth > _max_depth)
    _max_depth = tree->_traversed_depth;
}

void Range_Tree::add_sphere_to_node( One_Tree *tree, Tree_Node *&node, size_t isphere )
{
  ++tree->_current_level;
  if ( tree->_traversed_depth < tree->_current_level )
    tree->_traversed_depth = tree->_current_level;
  
  // do we need a new node?
  if (!node)
  {
    node = new_tree_node(tree);
    node->_sphere = isphere;
    // avoid creating so many size == 1 subtrees, delay creating a subtree until we are not a leaf
    //    if (tree->_d + 1 < num_dim())
    //    {
    //      node->_subtree = new_one_tree( tree->_d + 1 );
    //    }
  }
  // branch left or right
  else
  {
    const double &node_coord   = _spheres[ node->_sphere ][tree->_d];
    const double &sphere_coord = _spheres[ isphere       ][tree->_d];

    // if this was a leaf, create a new subtree and add self to it
    if (node->is_leaf() && tree->_d + 1 < num_dim())
    {
      node->_subtree = new_one_tree( tree->_d + 1 );
      add_sphere_to_one_tree( node->_subtree, node->_sphere);
    }

    // left/right
    if ( sphere_coord < node_coord )
    {
      add_sphere_to_node( tree, node->_left, isphere );
      // if this used to be a leaf, then add self to right
      if ( !node->_right )
      {
        add_sphere_to_node( tree, node->_right, node->_sphere);
        // my biggest item to the left is the new sphere
        node->_sphere = isphere;
        
        // new subtree, add both spheres
      }
    }
    else
    {
      add_sphere_to_node( tree, node->_right, isphere );
      // if this used to be a leaf, then add self to left
      if ( !node->_left )
      {
        add_sphere_to_node( tree, node->_left, node->_sphere);
        // my bigget item to the left is still myself
        // node->_sphere = unchanged
      }
    }
  }

  // recurse by d, on subtree
  if (tree->_d + 1 < num_dim())
    add_sphere_to_one_tree( node->_subtree, isphere);
  
  --tree->_current_level;
}

bool Range_Tree::needs_rebalance() const
{
  if (size() < 16)
    return false;

  // how far out of balance are we?
  // perfect is log2( size )
  // ok is log size * log log size
  size_t sz = size();
  size_t mylog(0);
  while ( sz >>= 1 ) ++mylog;
  size_t m( mylog );
  size_t myloglog(0);
  while ( m >>= 1 ) ++myloglog;
  
  const size_t perfect_depth = mylog + 1; // e.g. #nodes = 2^k - 1, then log is k-1
  const size_t ok_depth = perfect_depth * (size_t) ceil( sqrt(perfect_depth) ); // or * myloglog;
  // to do: also put a cap on how often rebalancing is done in the calling loop
  // debug
  if (_max_depth > ok_depth)
    std::cout << "Rebalancing depth:" << _max_depth << " size:" << size() << " perfect_depth: " << perfect_depth << " OK_depth:" << ok_depth << std::endl;
  return (_max_depth > ok_depth);
}


void Range_Tree::print( std::string name, std::ostream &out )
{
  std::cout << "======== Printing Range_Tree " << name << " ========" << std::endl;
  Sphere_Container::print( name, out );
#ifndef DEBUG_RANGE_TREE
  std::cout << "Tree information is limited. Re-compile with DEBUG_RANGE_TREE defined in Range_Tree.hpp" << std::endl;
#endif
  std::cout << std::endl << "Root" << std::endl;
  /* size_t nn = */ print_one_tree( &_tree );
  std::cout << "======== Done Printing Range_Tree " << name << " ========" << std::endl << std::endl;
}

size_t Range_Tree::print_one_tree( One_Tree *tree )
{
  std::cout << std::endl << "==== TREE ====" << std::endl;
#ifdef DEBUG_RANGE_TREE
  tree->_print_traversed_depth = 0;
  tree->_print_current_level = 0;
#endif
  size_t nn = print_node( tree->_root, tree );
  std::cout << "Tree had " << nn << " nodes (" << ((nn+1)/2) << " leaves)";
#ifdef DEBUG_RANGE_TREE
  assert( nn == tree->_max_node_id );
#endif
  if (tree == &_tree)
  {
    std::cout << " for " << size() << " spheres";
    //assert( num_nodes() == sum of nn );
    assert( nn + 1 == 2 * size() );
  }
#ifdef DEBUG_RANGE_TREE
  std::cout << ", and " << tree->_print_traversed_depth - 1 << " depth (starting at 0)." << std::endl;
#else
  std::cout << std::endl;
#endif

  std::cout << std::endl << "==== SUBTREES ====" << std::endl;
  print_subtrees( tree->_root, tree );
  return nn;
}

void Range_Tree::print_single_node( const Tree_Node *n, const One_Tree *tree ) const
{
  if (n)
  {
#ifndef DEBUG_RANGE_TREE
    std::cout << "Node:" << get_id(n) << " level:" << "?";
#else
    std::cout << "Node:" << get_id(n) << " level:" << n->_level;
#endif
    std::cout << " sphere:" << n->_sphere << " coord:";
    if ( tree )
      std::cout << spheres()[ n->_sphere ][tree->_d];
    else
      std::cout << " unknown";
  }
  else
  {
    std::cout << "Node: NULL ";
  }
  
  if (tree)
  {
#ifndef DEBUG_RANGE_TREE
    std::cout << " tree:" << tree;
#else
    std::cout << " tree:" << tree->_id;
#endif
    std::cout << " d:" << tree->_d;
  }
  else
  {
    std::cout << " tree: NULL";
  }
}

size_t Range_Tree::print_node(const Tree_Node *n, One_Tree *tree )
{
#ifdef DEBUG_RANGE_TREE
  if (++tree->_print_current_level > tree->_print_traversed_depth)
    tree->_print_traversed_depth = tree->_print_current_level;
#endif
  using std::cout;
  print_single_node( n, tree );
  size_t count(1);
  if (!n->_right)
    cout << "  LEAF. " << std::endl;
  else
  {
    std::cout << "  left ->" << get_id( n->_left ) << " and right->" << get_id( n->_right) << "." << std::endl;
    count += print_node( n->_left, tree );
    count += print_node( n->_right, tree );
  }
#ifdef DEBUG_RANGE_TREE
  --tree->_print_current_level;
#endif
  return count;
}

void Range_Tree::print_subtrees(const Tree_Node *n, One_Tree *tree )
{
  using std::cout;
  print_single_node( n, tree );
  if (!n->_subtree)
    cout << " no subtree. " << std::endl;
  else
  {
    cout << " subtree-> " << get_id(n) << std::endl;
#ifdef DEBUG_RANGE_TREE
    size_t nn =
#endif
    print_one_tree( n->_subtree );
#ifdef DEBUG_RANGE_TREE
    if ( n == tree->_root )
    {
      // assert( tree->_max_node_id == tree->num_nodes() );
      assert( nn == tree->_max_node_id );
    }
#endif
  }
  if (n->_left)
  {
    print_subtrees( n->_left, tree );
    assert( n->_right );
    print_subtrees( n->_right, tree );
  }
}

