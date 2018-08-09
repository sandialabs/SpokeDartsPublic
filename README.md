# SpokeDarts
SpokeDarts sphere-packing sampling in any dimension. Advancing front sampling from radial lines (spokes) through prior samples.
Two-spokes provides superior blue noise to traditional Maximal Poisson-disk Sampling, avoiding the distribution spike at exactly the sampling radius. Two-spokes is based on Two-radii from Section 3 in the paper "Variable Radii Poisson-Disk Sampling."
http://www.cs.sandia.gov/~samitch/bibliography_2007.html#var-radius-CCCG2012


## Compiling
The software has been compiled using Xcode on Mac OS X, and cmake and gcc on Linux. The code uses c++11 extensions. If you want to compile it multithreaded, you should consider disabling the use of the Mem_Pool memory managed pool of objects. See source/build/Readme.txt

## Interface
There is a basic interface taking input from standard input. See source/build/Readme.txt for instructions and infile.txt for an example. You can also edit the source code. For examples of how to run the various sampling algorithms and their options, see main.cpp. 

## Copyright and Credit

Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
SCR#:2084.0

To credit this work, please cite the following paper and software:


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

@misc{spokedartspubliccode,
 author={Muhammad A Awad and Mohamed S Ebeida and Scott A. Mitchell and  Anjul Patney and Ahmad A Rushdi and Laura P Swiler},
 title={{SpokeDartsPublic} Open-source Software},
 howpublished={v. 1.0, \url{https://github.com/samitch/SpokeDartsPublic}},
 year={2016}
}

% use the  "archivePrefix", "eprint", and "primaryClass" fields if your bibliography style handles it, 
% otherwise usepackage{url} in the latex file and use the "note" field in the bib file entry
% see also http://arxiv.org/hypertex/bibstyles/
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
% date = 5 aug 2014

