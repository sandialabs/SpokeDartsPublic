# build spokedarts using cmake (or xcode on mac)
# cd build (this directory)
cmake ..
make -j16

# At this point you should see the executable file spokedarts in this directory

# try to make a simple sample
./spokedarts < infile.txt
# 2
# 0
# 0.2
# 1
# 1
# 1
# 1
# 1
# 0
# 0

# sending infile.txt to standard input answers the questions the following way: 

Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
SCR#:2084.0

The domain is a unit cube
Enter the domain dimension
2
Good value entered: 2
Enter the domain periodicity (0=false, 1=true)
0
Good value entered: 0
Enter the sampling radius, i.e. the minimum intersample distance, between 0 and the box diagonal 1.414214
0.2
Good value entered: 0.2

Enter which algorithms you want to run
Bridson annular point sampling? (0/1) 
1
Good value entered: 1
Line-spokes? (0/1) 
1
Good value entered: 1
Favored-spokes? (0/1) 
1
Good value entered: 1
Two-spokes? (0/1) 
1
Good value entered: 1
Use the default random seed of 293879582987234? (0/1)
1
Good value entered: 1

Specify what kind of data structure you want to use to search for nearby samples.
The implemented options are a k-d tree, or exhaustive search.
Which is best depends on the number of samples and the dimension.
For middling dimensions and many points, a k-d tree is best.
For high dimensions and few points, exaustive search is best.
Because the domain dimension is 2, the default is a k-d tree
Enter the type of search. 0=keep default, 1=array, 3=k-d tree
0
Good value entered: 0
OK, we'll use a k-d tree
Do you want all searches to be global? Most users will say no=0. (0/1)
0
Good value entered: 0

========

You should see a bunch of files in this directory now, with names like
2_n_5_0spoke_0_mps_spheres.dat
2_n_5_0spoke_0_performance_final.txt
2_n_5_1favored_80_mps_spheres.dat
2_n_5_1favored_80_performance_final.txt
2_n_5_2two_112_mps_spheres.dat
2_n_5_2two_112_performance_final.txt
2_n_5_3Bridson_796_mps_spheres.dat
2_n_5_3Bridson_796_performance_final.txt

The leading "2" means the files are 2-dimensional
0spoke means the algorithm use was line-spokes
1favored is favored-spokes
2two is two-spokes
3Bridson is annular point sampling

The mps_spheres.dat files are text files of the following format
dimension  number_of_sample_spheres
bounding_box_min_coordinate[0] bounding_box_min_coordinate[1] ... bounding_box_min_coordinate[dimension]
bounding_box_max_coordinate[0] bounding_box_max_coordinate[1] ... bounding_box_max_coordinate[dimension]
1st_sample_coordinate[0] 1st_sample_coordinate[1] ... 1st_sample_coordinate[dimension] 1st_sample_radius
2nd_sample_coordinate[0] 1st_sample_coordinate[1] ... 2nd_sample_coordinate[dimension] 2nd_sample_radius
.
.
.
nth_sample_coordinate[0] nth_sample_coordinate[1] ... nth_sample_coordinate[dimension] 2nd_sample_radius

If you look in the source code, you will find options for outputing additional files and formats.
Darts_IO.hpp and calls from MPS.cpp, MPS::output_results, depending on the values of the _control_var values.


 

