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
#ifndef GRAMMATRIX_H
#define GRAMMATRIX_H
#include <vector>
#include <iostream>
#include <fstream>
#include <cfloat>
#include "matrix.h"
#include "DTW.h"
#include "KDTW.h"

/*
 * Gram matrix functions
 */
using namespace std;

double **vecmat2doubleptr(matrix<double> M);

matrix<double>  ComputeDTWGramMatrix(std::vector<UniformTimeseries> ds);
matrix<double>  ComputeDTWGramMatrix(std::vector<UniformTimeseries> ds, int nclass);

matrix<double>  ComputeKDTWGramMatrix(std::vector<UniformTimeseries> ds, double sigma);
matrix<double>  ComputeKDTWGramMatrix(std::vector<UniformTimeseries> ds, double sigma, int nclass);
matrix<double> ComputeKDistanceMatrix(std::vector<UniformTimeseries> ds, double sigma, int nclass);

#endif
