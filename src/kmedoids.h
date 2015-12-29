/*************************************************************************
eKATS0.1 (source code generated 2015-12-15)
Copyright (c) Pierre-François Marteau (eKATS project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
#ifndef KMEDOIDS_H
#define KMEDOIDS_H

#include <iostream>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include "UniformTimeseries.h"
#include "DTW.h"
#include "GramMatrix.h"
#include "cluster.h"

int getMedoid_min(matrix<double>  distMat, int *clusterid, int nc);
int getMedoid_max(matrix<double>  GramMat, int *clusterid, int nc);
std::vector<UniformTimeseries>  kmedoids_dtw(std::vector<UniformTimeseries> dsTrain, int nclass, int nclust);
std::vector<UniformTimeseries>  kmedoids_kdtw(std::vector<UniformTimeseries> dsTrain, int nclass, int nclust, double sigma);

#endif
