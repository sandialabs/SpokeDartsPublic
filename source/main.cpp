// main.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// SpokeDarts v. 1.0 software authors are Muhammad A. Awad, Mohamed S. Ebeida, Scott A. Mitchell, Ahmad A. Rushdi, and Laura P. Swiler.
// Additional contributions by Anjul Patney.

// To credit this work, please cite the following paper:

/*
@article{spokedarts_journal,
  author = {Scott A. Mitchell and Mohamed S. Ebeida and Muhammad A. Awad and Chonhyon Park and Anjul Patney and Laura P. Swiler and Dinesh Manocha and Li-Yi Wei and Ahmad A. Rushdi},
  title = {Efficient Blue Noise Sampling in any Dimension},
  journal = {submitted to ACM Trans. Graph.},
  volume = {},
  number = {},
  month = {},
  year = {2016},
  issn = {},
  pages = {},
  articleno = {},
  numpages = {},
  url = {},
  doi = {},
  acmid = {},
  publisher = {ACM},
  address = {New York, NY, USA},
  keywords = {Sampling; Algorithms, Theory, Experimentation; line sampling, high dimension, step blue noise, Delaunay graph, global optimization, motion planning},
  note={Open source software available from \url{https://github.com/samitch/SpokeDartsPublic}  An earlier version of this paper appears as \url{http://arxiv.org/abs/1408.1118}}
}
*/

/*
 % use the  "archivePrefix", "eprint", and "primaryClass" fields if your bibliography style handles it,
 % otherwise usepackage{url} in the latex file and use the "note" field in the bib file entry
 % see also http://arxiv.org/hypertex/bibstyles/
 @Article{spokedarts_arxiv,
	author = {Mohamed S. Ebeida and Scott A. Mitchell and Muhammad A. Awad and Chonhyon Park and Laura P. Swiler and Dinesh Manocha and Li-Yi Wei},
	title = {Spoke Darts for Efficient High Dimensional Blue Noise Sampling}
 journal   = "pre-print",
 volume    = "",
 year      = "2014",
 pages     = "12",
 
 eprinttype={arxiv},
 archivePrefix = "arXiv",
 eprint        = "1408.1118",
 primaryClass  = "cs.GR",
 version = {1},
 date={2014-8-5},
 note = {arXiv:1408.1118 [cs.GR] \url{http://arxiv.org/abs/1408.1118}},
	_note={arxiv submission 1038046},
 }
 % date = 5 aug 2014
*/


///////////////

#include "MPS_Spoke.hpp"
#include "MPS_Favored.hpp"
#include "MPS_Two_Spoke.hpp"
#include "MPS_Bridson.hpp"

#include "Scaling_Study.hpp"

int main(int argc, char *argv[])
{
  std::cout << "Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.\nSCR#:2084.0" << std::endl;
  
  // run the four sampling algorithms on the same domain
  if (/* DISABLES CODE */ (1))
  {
    
    // ==== domain dimension
    size_t num_dim(4);
    
    // ==== domain periodicity
    bool is_periodic = false;
    
    // ==== radius of spheres in packing
    double r(0.1);
    
    // radius limits for periodic domains, using implicit ghosts and full searches
    // there are other options where the only limit is r < 1
    // required: r < 0.161    for    line
    // required: r < 0.108    for favored
    // required: r < 0.071    for     two
    //
    
    // commonly used values
    //  num_dim = 2;
    //  r = 0.01;
    //  r = 0.03;
    //  r = 0.005;
    //  r = 1e-3;
    
    //  num_dim = 4;
    //  r = 0.1;
    //  r = 0.08;
    //  r = 0.05;
    
    //  num_dim = 5;
    //  r = 0.1;
    //  r = 0.16;
    //  r = 0.04;
    //  r = 0.02;
    
    //  num_dim = 6;
    //  r = 0.16;
    
    //  num_dim = 10;
    //  r = 0.37; // 0.4; // 0.1; // 0.2;
    
    
    // ==== random number generator
    
    // use your own rng
    Random my_random(num_dim);
    my_random.initiate_random_generator(293879582987234);
    
    // use the global static one
    // Random &my_random = Random::random_instance;
    
    //====== Bridson's annular point sampling
    
    MPS_Bridson mps_bridson;
    mps_bridson.create(num_dim, is_periodic, r, Search_Factory::TREE, false, false);
    
    //====== line-spokes
    
    MPS_Spoke mps_spoke;
    mps_spoke.set_rng(&my_random);
    
    if (/* DISABLES CODE */ (1))
    {
      // create the sampling using defaults
      mps_spoke.create(num_dim, is_periodic, r);
    }
    else
    {
      // fine-tune internal algorithms

      // what search routines should be used for proximity checking?
      // e.g. grids, array=exhaustive search, k-d trees, range trees
      // default is k-d tree
      // range trees may have bugs
      
      Search_Factory::Search_Type search_type = Search_Factory::UNSPECIFIED;
      bool all_searches_global = false;
      bool perm_ghosts = false; // determines how searches beyond periodic boundaries are handled
      
      // for high dimensions, >9, use an array and global searches always
      // Search_Factory::Search_Type search_type = Search_Factory::ARRAY;
      // bool all_searches_global = true;
      // bool perm_ghosts = false; // determines how searches beyond periodic boundaries are handled
    
      // if true, just create points in a small neighborhood around the central point,
      // for efficient testing in high dimensions of the value of beta that is achieved.
      mps_spoke._do_central_beta = false;

      // create the sampling
      mps_spoke.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
    }
    
    
    //====== favored spokes
    
    // examples
    // num_dim = 2; r = 0.007;
    // num_dim = 3; r = 0.015;
    // num_dim = 3; r = 0.02;
    
    MPS_Favored mps_favored;
    mps_favored.create(num_dim, is_periodic, r, Search_Factory::TREE);
    
    //===== two spokes
    
    // examples
    // num_dim = 2; r = 0.005;
    // num_dim = 3; r = 0.01;
    // num_dim = 4; r = 0.05;
    //  num_dim = 5; r = 0.071;
    
    MPS_Two_Spoke mps_two_spoke;
    mps_two_spoke.create(num_dim, is_periodic, r);
  }

  
  // study the value of achieved beta
  if (/* DISABLES CODE */ (0))
  {
    //    beta_study(6, true, 4);
    beta_study2(6, true);
  }
  // study runtime and memory and distribution size perforfance across r, dimension, and algorithm
  if (/* DISABLES CODE */ (0))
  {
    //    scaling_study_1();
    //    scaling_study_2();
    //    scaling_study_3();
    //    scaling_study_4();
    //    scaling_study_5();
    //    scaling_study_6();
    //    scaling_study_7();
    //    scaling_study_8();
    //    scaling_study_9();
    //    scaling_study_10();
    //    scaling_study_11();
    //    scaling_study_12();
    scaling_study_bridson_1();
  }
  
  // no errors
  return 0;
}
