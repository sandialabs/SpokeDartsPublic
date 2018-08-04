//  Search_Factory.cpp
//  spokes
//

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0


#include "Search_Factory.hpp"
#include "Search_Array.hpp"
#include "Search_Ghost_Array.hpp"
#include "Search_Grid.hpp"
#include "Search_Tree.hpp"
#include "Search_Range.hpp"


const char* Search_Factory::search_description[] = { "unspecified", "exhaustive search array, O(n)", "a uniform background grid", "a k-d tree", "a range tree (experimental!)" };
const char* Search_Factory::search_name[] = { "unspecified", "array", "grid", "k-d tree", "range tree (experimental!)" };


Search_Structure* Search_Factory::new_search( Search_Type search_type, Spheres *spheres, bool is_global, double search_distance, double xmax, double xmin )
{
  Search_Structure *search(0);

  switch (search_type)
  {
    case GRID:
    {
      search = new Search_Grid(*spheres, search_distance, xmax, xmin);
    }
    break;

    case RANGE:
    {
      search = new Search_Range(*spheres);
    }
    break;
      
    case TREE:
    case UNSPECIFIED:
    {
      search = new Search_Tree(*spheres);
    }
    break;

    case ARRAY:
    default:
    {
      Search_Array *search_array = new Search_Array(*spheres, is_global);
      search = search_array;
      
      // check to see if we should use a global ghost array instead
      Ghost_Spheres *ghosts = spheres->cast_to_ghost();
      if ( is_global && ghosts )
        search = new Search_Ghost_Array( *ghosts, search_array );
    }
    break;
      
  }
  
  return search;
}
