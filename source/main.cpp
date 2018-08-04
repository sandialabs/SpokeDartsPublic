// main.cpp

// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0

// SpokeDarts v. 1.0 software authors are Muhammad A. Awad, Mohamed S. Ebeida, Scott A. Mitchell, Ahmad A. Rushdi, and Laura P. Swiler.
// Additional contributions by Anjul Patney.

// To credit this work, please cite the following paper:
/*
 @article{Mitchell:2018:SHB:3191713.3194657,
 author = {Scott A. Mitchell and Mohamed S. Ebeida and Muhammad A. Awad and Chonhyon Park and Anjul Patney and Ahmad A. Rushdi and Laura P. Swiler and Dinesh Manocha and Li-Yi Wei},
 title = {Spoke-Darts for High-Dimensional Blue-Noise Sampling},
 journal = {ACM Trans. Graph.},
 volume = {37},
 number = {2},
 month = {May},
 year = {2018},
 issn = {0730-0301},
 pages = {22:1--22:20},
 articleno = {22},
 numpages = {20},
 url = {http://doi.acm.org/10.1145/3194657},
 doi = {10.1145/3194657},
 acmid = {3194657},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {Delaunay graph, Line sampling, blue noise, global optimization, high dimension, motion planning},
 note={Open source software available from \url{https://github.com/samitch/SpokeDartsPublic}}
 }
 % may 2018
 */

// to cite the software itself
/*
 @misc{spokedartspubliccode,
 author={Muhammad A Awad and Mohamed S Ebeida and Scott A. Mitchell and  Anjul Patney and Ahmad A Rushdi and Laura P Swiler},
 title={{SpokeDartsPublic} Open-source Software},
 howpublished={v. 1.0, \url{https://github.com/samitch/SpokeDartsPublic}},
 year={2016}
 }
*/

// this is an arxiv version of the above ACM Trans. Graph. paper
/*
 % use the  "archivePrefix", "eprint", and "primaryClass" fields if your bibliography style handles it,
 % otherwise usepackage{url} in the latex file and use the "note" field in the bib file entry
 % see also http://arxiv.org/hypertex/bibstyles/
 % The paper password for 1408.1118 is sgwsf
 @Article{spokedarts_arxiv,
 author = {Scott A. Mitchell and Mohamed S. Ebeida and Muhammad A. Awad and Chonhyon Park and Anjul Patney and Ahmad A. Rushdi and Laura P. Swiler and Dinesh Manocha and Li-Yi Wei},
 title = {Spoke-Darts for High-Dimensional Blue-Noise Sampling}
 journal   = "pre-print",
 volume    = "",
 year      = "2018",
 pages     = "19",
 
 eprinttype={arxiv},
 archivePrefix = "arXiv",
 eprint        = "1408.1118",
 primaryClass  = "cs.GR",
 version = {3},
 date={2018-6-13},
 note = {arXiv:1408.1118 [cs.GR] \url{http://arxiv.org/abs/1408.1118}},
 _note={arxiv submission 1038046},
 }
*/


///////////////

#include "MPS_Spoke.hpp"
#include "MPS_Favored.hpp"
#include "MPS_Two_Spoke.hpp"
#include "MPS_Bridson.hpp"

#include "Scaling_Study.hpp"

// read a number from standard input
template <class Num>
void ui_num(std::string s, Num &n)
{
  bool good_value = false;
  do
  {
    std::cout << s << std::endl;
    std::cin >> n;
    good_value = std::cin && !std::cin.fail() && !std::isnan(n);
    if (!good_value)
    {
      std::cin.clear();
      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      std::cout << "Bad value entered: " << n << std::endl;
    }
  } while (!good_value);
  std::cout << "Good value entered: " << n << std::endl;
}

// get domain and sampling parameters from crude user interface
void ui_parameters(size_t &num_dim, bool &is_periodic, double &r)
{
  std::cout << "\nThe domain is a unit cube" << std::endl;
  do
    ui_num("Enter the domain dimension", num_dim);
  while (num_dim==0 && num_dim>10000);
  ui_num("Enter the domain periodicity (0=false, 1=true)", is_periodic);
  const double diagonal_length = sqrt(num_dim);
  do
    ui_num("Enter the sampling radius, i.e. the minimum intersample distance, between 0 and the box diagonal " + std::to_string(diagonal_length), r);
  while ( r <= 0. && r > diagonal_length );
}

int main(int argc, char *argv[])
{
  std::cout << "Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.\nSCR#:2084.0" << std::endl;
  
  if (/* DISABLES CODE */ (1))
  {
    size_t num_dim(4);
    bool is_periodic = false;
    double r(0.1);
    ui_parameters(num_dim,is_periodic,r);
    
    std::cout << "\n" << "Enter which algorithms you want to run" << std::endl;
    bool do_bridson(false), do_line, do_favored, do_two;
    ui_num("Bridson annular point sampling? (0/1) ", do_bridson);
    ui_num("Line-spokes? (0/1) ", do_line);
    ui_num("Favored-spokes? (0/1) ", do_favored);
    ui_num("Two-spokes? (0/1) ", do_two);

    // change rng seed? or use default
    const unsigned long default_seed = 293879582987234;
    auto seed = default_seed;
    bool keep_seed;
    ui_num("Use the default random seed of " + std::to_string(default_seed) + "? (0/1)", keep_seed);
    if (!keep_seed)
      ui_num("Enter the random number generator seed", seed);
    
    // change search type?
    Search_Factory::Search_Type search_type = ((num_dim < 7) ? Search_Factory::TREE : Search_Factory::ARRAY);
    // if (num_dim < 4)
    //  search_type = Search_Factory::GRID;
    std::cout << "\nSpecify what kind of data structure you want to use to search for nearby samples.\n";
    std::cout << "The implemented options are a k-d tree, or exhaustive search.\n";
    //    std::cout << "(A uniform grid and a range tree are only partially implemented)\n";
    std::cout << "Which is best depends on the number of samples and the dimension.\n";
    //     std::cout << "For very low dimensions, a background grid is best.\n";
    std::cout << "For middling dimensions and many points, a k-d tree is best.\n";
    std::cout << "For high dimensions and few points, exaustive search is best.\n";
    Search_Factory sf;
    std::cout << "Because the domain dimension is " << num_dim << ", the default is " << sf.search_description[ (size_t) search_type] << std::endl;
    int search_id;
    do
    {
      ui_num("Enter the type of search. 0=keep default, 1=" + std::string(sf.search_name[1]) + ", 3=" + std::string(sf.search_name[3]), search_id);
    } while (search_id < 0 || search_id > 3 || search_id == 2);
    if (search_id>0)
      search_type = (Search_Factory::Search_Type) search_id;
    std::cout << "OK, we'll use " << sf.search_description[ (size_t) search_type] << std::endl;
    
    bool all_searches_global;
    ui_num("Do you want all searches to be global? Most users will say no=0. (0/1)", all_searches_global);
    bool perm_ghosts = false;

    if (do_bridson)
    {
      MPS_Bridson mps;
      Random my_random(num_dim);
      my_random.initiate_random_generator(seed);
      mps.set_rng(&my_random);
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
    }

    if (do_line)
    {
      MPS_Spoke mps;
      Random my_random(num_dim);
      my_random.initiate_random_generator(seed);
      mps.set_rng(&my_random);
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
    }

    if (do_favored)
    {
      MPS_Favored mps;
      Random my_random(num_dim);
      my_random.initiate_random_generator(seed);
      mps.set_rng(&my_random);
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
    }

    if (do_two)
    {
      MPS_Two_Spoke mps;
      Random my_random(num_dim);
      my_random.initiate_random_generator(seed);
      mps.set_rng(&my_random);
      mps.create(num_dim, is_periodic, r, search_type, all_searches_global, perm_ghosts);
    }
    
  }
  
  // here are examples of hard-coded runs, that you can modify for specific studies
  // run the four sampling algorithms on the same domain
  if (/* DISABLES CODE */ (0))
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
