//  ScalingStudy.h
//  spokes
//
// Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
// SCR#:2084.0
//

#ifndef __spokes__ScalingStudy__
#define __spokes__ScalingStudy__

#include <stdio.h>

void scaling_study_1();
void scaling_study_2();
void scaling_study_3();
void scaling_study_4();
void scaling_study_5();

void scaling_study_6(); // extreme dimensions

void scaling_study_7();

void scaling_study_8(); // periodic array

void scaling_study_9(); // periodic array, using the new marching algorithm for trimming
void scaling_study_10(); // periodic array, using the new marching algorithm for trimming
void scaling_study_11(); // periodic array, using the new marching algorithm for trimming
void scaling_study_12(); // periodic array, using the new marching algorithm for trimming

void scaling_study_bridson_1();

// lots of trials, one point
void beta_study(size_t num_dim = 2, bool estimate_beta = false, size_t layer_limit = (size_t) -1);

// one trial, 10,000 points, estimate beta

// central beta, r is small and aperiodic, go out some distances
void beta_study2(size_t num_dim = 2, bool estimate_beta = true );

// central beta, but r is large and domain is periodic
void beta_study2( size_t num_dim, bool estimate_beta, size_t layer_limit );


#endif /* defined(__spokes__ScalingStudy__) */
