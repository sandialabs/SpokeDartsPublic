// Mem_Pool.hpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// block memory allocate generic class objects, and re-use individual items

#ifndef MEM_POOL_HPP
#define MEM_POOL_HPP

//=================================================
#ifdef MPS_MEMORY_MANAGEMENT

template< class Object_Type>
class Mem_Pool
{
  class Mem_Block
  {
  public:
    Mem_Block( size_t mem_size )
    {
      // avoid calling the constructor
      // _mem = new Object_Type[mem_size]; conceptually, but without the constructor
      void * raw_mem =  malloc( mem_size * sizeof(Object_Type) );
      _mem = (Object_Type*) raw_mem;
      _next_free = _mem;
      _last_mem = &_mem[mem_size]; // this memory location is unallocated, bad
    }
    ~Mem_Block() { delete [] _mem; }
    Object_Type *get_next_mem()
    {
      Object_Type *current_free = _next_free;
      ++_next_free;
      // if _next_free == _last_mem, then that is OK since the last byte of current_free was OK
      return (_next_free > _last_mem) ? 0 : current_free;
    }
  private:
    Object_Type * _mem;
    Object_Type * _next_free;
    Object_Type * _last_mem;
  };
  
  size_t _block_size;
  std::vector<Mem_Block*> _mem_blocks;
  
  // always succeeds
  Object_Type *get_block_mem()
  {
    Object_Type *mem(0);
    if (!_mem_blocks.empty())
      mem = _mem_blocks.back()->get_next_mem();
    if (!mem)
    {
      _block_size *= 2; // or 1.5?
      _mem_blocks.push_back( new Mem_Block( _block_size ) );
      mem = _mem_blocks.back()->get_next_mem();
    }
    assert(mem);
    return mem;
  }

  
  // pool of freed memory
  typedef std::vector<Object_Type *> Object_Pool;
  Object_Pool _object_pool;

  // return memory from the pool, if any, or 0
  Object_Type *get_pool_mem()
  {
    if (!_object_pool.empty())
    {
      Object_Type *reused_object = _object_pool.back();
      _object_pool.pop_back();
      return reused_object;
    }
    return 0;
  }

  
public:
  
  Mem_Pool( size_t block_size ) : _block_size( block_size ) {}
  ~Mem_Pool()
  {
    // ignore objects in the object pool; they are in some block
    _object_pool.clear();

    // free the blocks
    while (!_mem_blocks.empty())
    {
      delete _mem_blocks.back();
      _mem_blocks.pop_back();
    }
  }
  
  // always succeeds; the managed class's new should call this to get a new object
  Object_Type *get_mem()
  {
    Object_Type *mem = get_pool_mem();
    if (!mem)
      mem = get_block_mem();
    assert(mem);
    return mem;
  }

  // the managed class's new should call this when done
  void put_mem( void *p )
  {
    assert(p);
    _object_pool.push_back( (Object_Type*) p);
  }
  
};

/* For any class you want to memory manage

// idiom
// set up the static memory manager for this class
// in hpp
#include "Mem_Pool.hpp"
class Object_Type
{
  // static so we always use this pool
  static Mem_Pool<Object_Type> _mem_pool;

  // idiom
  void * operator new( size_t size )
  {
    assert( size == sizeof( Object_Type ));
    return _mem_pool.get_mem();
    // constructor will initialize values if desired
  }

  // idiom
  void operator delete( void *p )
  {
    if (p)
      _mem_pool.put_mem(p);
    // caller should set p = 0;
  }

};

// in cpp file
Object_Type::Mem_Pool<Object_Type> _mem_pool( 1000 ); // block size

*/

// end MPS_MEMORY_MANAGEMENT
#endif
//=================================================


// end _HPP
#endif
