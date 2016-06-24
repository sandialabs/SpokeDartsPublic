// Range_Tree.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// This is the data for the tree

#ifndef RANGE_TREE_HPP
#define RANGE_TREE_HPP

#include "Sphere_Container.hpp"
#include "Point_Tool.hpp"
#include "Sphere_Container_Array.hpp"
//=================================================
#ifdef MPS_MEMORY_MANAGEMENT
#include "Mem_Pool.hpp"
#endif
//=================================================

// uncomment this to get a lot of debugging info printed out
// #define DEBUG_RANGE_TREE

class Range_Tree : public Sphere_Container
{

public:

  Range_Tree(Spheres &spheres)
  : Sphere_Container( spheres ), _array( spheres ), _tree( 0 ), _max_depth(0), _tree_ispheres(0), _tree_ispheres_size(0)
  { new_memory(); };

  virtual ~Range_Tree()
  { delete_memory(); }

  // Sphere_Container interface functions

  // get rid of contents of the current tree
  virtual void clear();

  // add one sphere to a tree, and the array
  // becomes slow if the tree is unbalanced
  virtual void add_sphere(size_t isphere);
  virtual void add_spheres(Sphere_Container *spheres);
  
  // after add_sphere, one may check whether the tree needs rebalancing
  bool needs_rebalance() const;

  // rebalance the tree for more uniform depth
  virtual void rebalance();

  // number of spheres in the tree
  virtual size_t size() const {return _array.size();}

  // number of tree nodes that are in this Range_Tree,
  // summed over all One_Trees.
  // This should be the same whether the tree is balanced or not
  virtual size_t num_nodes() const 
  { return 2 * size() * num_dim(); }

  // access the array form of the tree
  const Sphere_Container_Array &array() { return _array; }

  // tree[i]
  //
  // rvalue access only through operator []
  // size_t node = tree[9];
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

  // print out the tree structure, in text
  virtual void print( std::string name = std::string(), std::ostream &out = std::cout);
  
protected:
  // methods

  // add node to tree
  // the node has to be added to the array separately, e.g. by the search_range_near_spheres method
  void add_sphere_to_tree(const size_t isphere);
  
  // memory allocation
  void new_memory();
  void delete_memory();

protected:

  // data
  struct One_Tree;
  struct Tree_Node
  {
    // if leaf, this is the sphere it contains
    // if not a leaf, then this is the sphere with the largest left value
    size_t _sphere;

    // if not leaf
    One_Tree *_subtree; // of the next dimension
    Tree_Node *_left, *_right;

    bool is_leaf()
    {       
      assert( is_wellformed() );
      return (_right == 0);
    }
    bool is_wellformed(One_Tree *tree = 0, size_t num_dim = 0 )
    {
      assert( this );
      // valid sphere?
      // id? level?
      // either both sides are null (leaf), or both sides are non-null (non-leaf)
      assert( (_left && _right) || (!_left && !_right) );
      // if it isn't a leaf, then this is either the last dimension, or there is a non-null subtree
      if (_left && _right)
      {
        assert( !tree || ((tree->_d + 1 == num_dim) && !_subtree) || ((tree->_d < num_dim) && _subtree) );
        return( !tree || ((tree->_d + 1 == num_dim) && !_subtree) || ((tree->_d < num_dim) && _subtree) );
      }
      // its a leaf
      if (!_left && !_right )
      {
        assert( !_subtree );
        return( !_subtree );
      }
      return false;
    }
    
    Tree_Node(): _left(0), _right(0), _subtree(0), _sphere( Spheres::bad_sphere_index() )
#ifdef DEBUG_RANGE_TREE
    , _id( 0 ), _level (0)
#endif
    {}

    // debug
#ifdef DEBUG_RANGE_TREE
    size_t _id;
    size_t _level;
#endif
    
    
//=================================================
#ifdef MPS_MEMORY_MANAGEMENT

    // set up the static memory manager for this class
    // static so we always use this pool
    static Mem_Pool<Tree_Node> _mem_pool;

    void * operator new( size_t size )
    {
      assert( size == sizeof( Tree_Node ));
      return _mem_pool.get_mem();
      // constructor will initialize values if desired
    }
    
    void operator delete( void *p )
    {
      if (p)
        _mem_pool.put_mem(p);
      // caller should set p = 0;
    }
#endif
//=================================================
  }; // end class Tree_Node


  Tree_Node * new_tree_node( One_Tree *tree );
  // recursively deletes the dependent objects
  void delete_tree_node( Tree_Node *node );

  // tree on one dimension
  struct One_Tree
  {
    Tree_Node *_root;  

    size_t _d; // dimension, index of coordinates of spheres

    // max depth
    size_t _traversed_depth;
    // current depth
    size_t _current_level;

    One_Tree(size_t d) : _root(0), _d(d), _traversed_depth(0), _current_level(0)
#ifdef DEBUG_RANGE_TREE
    , _id( ++_max_id ), _max_node_id(0)
#endif
    {}
    
#ifdef DEBUG_RANGE_TREE
    static size_t _max_id;
    size_t _max_node_id;
    size_t _id;
    
    size_t _print_current_level;
    size_t _print_traversed_depth;

#endif

//=================================================
#ifdef MPS_MEMORY_MANAGEMENT
    
    // set up the static memory manager for this class
    // static so we always use this pool
    static Mem_Pool<One_Tree> _mem_pool;
    
    void * operator new( size_t size )
    {
      assert( size == sizeof( One_Tree ));
      return _mem_pool.get_mem();
      // constructor will initialize values if desired
    }
    
    void operator delete( void *p )
    {
      if (p)
        _mem_pool.put_mem(p);
      // caller should set p = 0;
    }
#endif
//=================================================
  }; // end class One_Tree
  

  One_Tree * new_one_tree( size_t d );
  // this takes care of recursively deleting the dependent objects
  void delete_one_tree( One_Tree *tree );

  One_Tree _tree; // top level tree, _d = 0

  bool last_dimension( One_Tree *tree )
  {
    assert( tree->_d < num_dim() );
    return tree->_d + 1 == num_dim();
  }

  // keep track of tree depth during insertions
  size_t _max_depth;
  
  // straight array of all the indices in the tree
  Sphere_Container_Array _array;

  // Speed up rebalance, by pre-allocating the arrays used for sorting
  size_t **_tree_ispheres;
  size_t _tree_ispheres_size;
  
  virtual void reset_iterator_it( size_t &it ) const
  {
    return _array.reset_iterator_it( it );
  }
  virtual size_t next_it( size_t &it ) const
  {
    return _array.next_it( it );
  }

private:

  void clear_tree();

  // recursive functions called by public non-recursive function
  void rebalance_tree(One_Tree *tree, const size_t *ispheres, size_t isize);
  void rebalance_node(One_Tree *tree, Tree_Node *&node, size_t *ispheres, size_t isize);

  void add_sphere_to_one_tree(One_Tree *tree, const size_t isphere);
  void add_sphere_to_node( One_Tree *tree, Tree_Node *&node, size_t isphere );

protected:
  // print
  // conceptually these are all const, but in debug mode tree has some level data that is set and printed
  // print nodes, then print subtrees of nodes. Return the nubmer of nodes printed (in the tree).
  size_t print_one_tree( One_Tree *tree );
  // dump node, then dump left then dump right. DFS
  size_t print_node( const Tree_Node *n, One_Tree *tree );
  // print subtree, then recurse on left and right
  void print_subtrees( const Tree_Node *n, One_Tree *tree );
  // output the actual text describing a node
  void print_single_node( const Tree_Node *n, const One_Tree *tree ) const;
  
#ifdef DEBUG_RANGE_TREE
  size_t get_id( const Tree_Node *n ) const
  { return n->_id; }
#else
  const void * get_id( const Tree_Node *n ) const
  { return n; }
#endif


};

#endif
